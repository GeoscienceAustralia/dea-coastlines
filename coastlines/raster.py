#!/usr/bin/env python
# coding: utf-8

# This code conducts raster generation for DE Africa Coastlines:

#     * Load stack of all available Landsat 5, 7 and 8 satellite imagery
#       for a location using ODC Virtual Products
#     * Convert each satellite image into a remote sensing water index
#       (MNDWI)
#     * For each satellite image, model ocean tides into a tide modelling
#       grid based on exact time of image acquisition
#     * Interpolate tide heights into spatial extent of image stack
#     * Mask out high and low tide pixels by removing all observations
#       acquired outside of 50 percent of the observed tidal range
#       centered over mean sea level
#     * Combine tidally-masked data into annual median composites from
#       representing the most representative position of the shoreline
#       at approximately mean sea level tide height each year (0 metres
#       Above Mean Sea Level).

import multiprocessing
import os
import sys
import warnings
from collections import Counter
from functools import partial

import click
import datacube
import geopandas as gpd
import numpy as np
import odc.algo
import otps
import pandas as pd
import pytz
import xarray as xr
from affine import Affine
from datacube.utils.cog import write_cog
from datacube.utils.geometry import CRS, GeoBox, Geometry
from datacube.utils.masking import make_mask
from datacube.virtual import catalog_from_file
from dea_tools.dask import create_local_dask_cluster
from dea_tools.spatial import interpolate_2d
from shapely.geometry import shape

from coastlines.utils import configure_logging, load_config

# Hide warnings
warnings.simplefilter(action="ignore", category=FutureWarning)


def hillshade(dem, elevation, azimuth, vert_exag=1, dx=30, dy=30):

    # For each supplied azimuth, compute hillshade
    from matplotlib.colors import LightSource

    hs = LightSource(azdeg=azimuth, altdeg=elevation).hillshade(
        dem, vert_exag=vert_exag, dx=dx, dy=dy
    )
    return hs


def terrain_shadow(x, elevation, threshold=0.4, radius=3):

    from skimage.morphology import binary_dilation, binary_opening, disk

    hs = hillshade(elevation, x.sun_elevation, x.sun_azimuth)
    hs = hs < threshold
    hs = binary_opening(hs, disk(radius))
    hs = binary_dilation(hs, disk(radius))

    return xr.DataArray(hs, dims=["y", "x"])


def sun_angles(dc, query):

    import numpy as np
    import datacube
    from datacube.api.query import query_group_by
    from datacube.model.utils import xr_apply

    # Customise query to match query used to load raster data
    query_subset = {k: v for k, v in query.items() if k not in ["dask_chunks"]}
    query_subset.update(
        product=["ls5_sr", "ls7_sr", "ls8_sr", "ls9_sr"],
        collection_category="T1",
        group_by="solar_day",
    )

    # Load the metadata from each product into an xarray after grouping by solar day
    gb = query_group_by(**query_subset)
    datasets = dc.find_datasets(**query_subset)
    dataset_array = dc.group_datasets(datasets, gb)
    sun_azimuth = xr_apply(
        dataset_array,
        lambda t, dd: np.mean([d.metadata.eo_sun_azimuth for d in dd]),
        dtype=float,
    )
    sun_elevation = xr_apply(
        dataset_array,
        lambda t, dd: np.mean([d.metadata.eo_sun_elevation for d in dd]),
        dtype=float,
    )

    # Combine into new xarray.Dataset
    sun_angles_ds = xr.merge(
        [sun_elevation.rename("sun_elevation"), sun_azimuth.rename("sun_azimuth")]
    )

    return sun_angles_ds


def mask_terrain_shadow(dc, query, ds):

    # Compute solar angles for all satellite image timesteps
    sun_angles_ds = sun_angles(dc, query)

    # Load DEM into satellite data geobox
    dem_ds = dc.load(product="dem_srtm", like=ds.geobox, resampling="cubic").squeeze(
        "time", drop=True
    )
    dem_ds = dem_ds.where(dem_ds.elevation >= 0)

    # Identify terrain shadow across all timesteps
    terrain_shadow_ds = multiprocess_apply(
        sun_angles_ds,
        dim="time",
        func=partial(terrain_shadow, elevation=dem_ds.elevation.values),
    )

    # Remove terrain shadow pixels from satellite data
    return ds.where(~terrain_shadow_ds)


def load_water_index(dc, query, yaml_path, product_name="ls_nbart_mndwi"):
    """
    This function uses virtual products to load Landsat 5, 7 and 8 data,
    calculate a custom remote sensing index, and return the data as a
    single xarray.Dataset.

    To minimise resampling effects and maintain the highest data
    fidelity required for subpixel coastline extraction, this workflow
    applies masking and index calculation at native resolution, and
    only re-projects to the most common CRS for the query using cubic
    resampling in the final step.

    Parameters:
    -----------
    dc : datacube.Datacube object
        Datacube instance used to load data.
    query : dict
        A dictionary containing query parameters passed to the
        datacube virtual product (e.g. same as provided to `dc.load`).
    yaml_path : string
        Path to YAML file containing virtual product recipe.
    product_name : string, optional
        Name of the virtual product to load from the YAML recipe.

    Returns:
    --------
    ds : xarray.Dataset
        An `xarray.Dataset` containing a time series of water index
        data (e.g. MNDWI) for the provided datacube query
    """
    import mock

    def custom_native_geobox(ds, measurements=None, basis=None):
        """
        Obtains native geobox info from dataset metadata
        """
        geotransform = ds.metadata_doc["grids"]["default"]["transform"]
        shape = ds.metadata_doc["grids"]["default"]["shape"]
        crs = CRS(ds.metadata_doc["crs"])
        affine = Affine(
            geotransform[0], 0.0, geotransform[2], 0.0, geotransform[4], geotransform[5]
        )

        return GeoBox(width=shape[1], height=shape[0], affine=affine, crs=crs)

    # Load in virtual product catalogue and select MNDWI product
    catalog = catalog_from_file(yaml_path)
    product = catalog[product_name]

    # Construct a new version of the product using most common CRS
    # Determine geobox with custom function to increase lazy loading
    # speed (will eventually be done automatically within virtual
    # products)
    with mock.patch(
        "datacube.virtual.impl.native_geobox", side_effect=custom_native_geobox
    ):

        # Identify most common CRS
        bag = product.query(dc, **query)
        crs_list = [str(i.crs) for i in bag.contained_datasets()]
        crs_counts = Counter(crs_list)
        crs = crs_counts.most_common(1)[0][0]

        # Pass CRS to product load
        settings = dict(
            output_crs=crs,
            resolution=(-30, 30),
            align=(15, 15),
            resampling={"pixel_quality": "nearest", "*": "cubic"},
        )
        box = product.group(bag, **settings, **query)
        ds = product.fetch(box, **settings, **query)

    # Rechunk if smallest chunk is less than 10
    if ((len(ds.x) % 3000) <= 10) or ((len(ds.y) % 3000) <= 10):
        ds = ds.chunk({"x": 3200, "y": 3200})

    # Identify pixels that are either cloud, cloud shadow or nodata
    nodata = make_mask(ds["pixel_quality"], nodata=True)
    mask = (
        make_mask(ds["pixel_quality"], cloud="high_confidence")
        | make_mask(ds["pixel_quality"], cloud_shadow="high_confidence")
        | nodata
    )

    # Apply opening to remove long narrow false positive clouds along
    # the coastline, then dilate to restore cloud edges
    mask_cleaned = odc.algo.mask_cleanup(
        mask, mask_filters=[("opening", 20), ("dilation", 5)]
    )
    ds = ds.where(~mask_cleaned & ~nodata)

    # Mask any invalid pixel values outside of 0 and 1
    green_bool = (ds.green >= 0) & (ds.green <= 1)
    swir_bool = (ds.swir_1 >= 0) & (ds.swir_1 <= 1)
    nir_bool = (ds.nir >= 0) & (ds.nir <= 1)
    ds = ds.where(green_bool & swir_bool)

    # Compute MNDWI
    ds[["mndwi"]] = (ds.green - ds.swir_1) / (ds.green + ds.swir_1)
    ds[["ndwi"]] = (ds.green - ds.nir) / (ds.green + ds.nir)

    # Apply terrain mask to remove deep shadows that can be
    # be mistaken for water
    ds = mask_terrain_shadow(dc, query, ds)

    return ds[["mndwi", "ndwi"]]


def otps_tides(lats, lons, times, timezone=None):
    """
    Model tide heights for one or more locations and times using the
    OTPS TPXO8 tidal model.

    Parameters:
    -----------
    lats, lons : numeric or list of numeric values
        One or more latitudes and longitude coordinates used to define
        the location at which to model tides.
    times : datetime.datetime or list of datetime.datetimes
        One or more `datatime.datetime` objects providing the times at
        which to model tides. By default these are assumed to be in UTC
        time; if this is not the case, use `timezone` below.
    timezone : string, optional
        If provided `datatime.datetime`s are not in UTC times, use this
        parameter to declare a timezone. E.g. to model tides for times
        expressed in local time at Darwin, Australia, provide
        `timezone='Australia/Darwin'`. Defaults to `None`, which assumes
        provided times are UTC. This is used to convert all times to UTC
        using the `pytz` module. For a full list of timezones, run:
        `import pytz; pytz.all_timezones`.

    Returns:
    --------
    tidepoints_df : pandas.DataFrame
        An `pandas.DataFrame` with a "time" index, "lat" and "lon"
        columns, and a "tide_m" column giving tide heights at each
        point location.
    """

    # Convert to list if provided as individual values
    if not isinstance(lats, list):
        lats = [lats]
    if not isinstance(lons, list):
        lons = [lons]
    if not isinstance(times, list):
        times = [times]

    # If a timezone is provided, localise the input times then
    # standardise to UTC times
    if timezone:
        times = [
            pytz.timezone(timezone).localize(time).astimezone(pytz.utc)
            for time in times
        ]

    # Create list of lat/lon/time scenarios to model tides
    observed_timepoints = [
        otps.TimePoint(lon, lat, time) for time in times for lon, lat in zip(lons, lats)
    ]

    # Model tides for each lat/lon/time
    observed_predictedtides = otps.predict_tide(observed_timepoints)

    # Output results into pandas.DataFrame
    tidepoints_df = pd.DataFrame(
        [
            (i.timepoint.timestamp, i.timepoint.lon, i.timepoint.lat, i.tide_m)
            for i in observed_predictedtides
        ],
        columns=["time", "lon", "lat", "tide_m"],
    )

    return tidepoints_df.set_index("time")


def model_tide_points(ds, points_gdf, extent_buffer=0.05):
    """
    Takes an xarray.Dataset (`ds`), extracts a subset of tide modelling
    points from a geopandas.GeoDataFrame based on`ds`'s extent, then
    uses the OTPS tidal model to model tide heights for every point
    at every time step in `ds`.

    The output is a geopandas.GeoDataFrame with a "time" index
    (matching the time steps in `ds`), and a "tide_m" column giving the
    tide heights at each point location.

    Parameters:
    -----------
    ds : xarray.Dataset
        An `xarray.Dataset` containing a time series of water index
        data (e.g. MNDWI) for the provided datacube query.
    points_gdf : geopandas.GeoDataFrame
        A `geopandas.GeoDataFrame` containing spatial points used to
        model tides using the OTPS tidal model.
    extent_buffer : float, optional
        A float giving the extent in degrees to buffer the satellite
        imagery dataset (`ds`) when selecting spatial points used
        to model tides. This buffer creates overlap between analysis
        areas, which ensures that modelled tides are seamless when
        clipped back to the dataset extent in a subsequent step.

    Returns:
    --------
    tidepoints_gdf : geopandas.GeoDataFrame
        An `geopandas.GeoDataFrame` containing modelled tide heights
        with an index based on each timestep in `ds`.
    """

    # Obtain extent of loaded data, and buffer to ensure that tides
    # are modelled reliably and comparably across grid tiles
    ds_extent = shape(ds.geobox.geographic_extent.json)
    buffered = ds_extent.buffer(extent_buffer)
    subset_gdf = points_gdf[points_gdf.geometry.intersects(buffered)]

    # Extract lon, lat from tides, and time from satellite data
    x_vals = subset_gdf.geometry.x.tolist()
    y_vals = subset_gdf.geometry.y.tolist()
    observed_datetimes = ds.time.data.astype("M8[s]").astype("O").tolist()

    # Model tides for each coordinate and time
    tidepoints_df = otps_tides(lats=y_vals, lons=x_vals, times=observed_datetimes)

    # Convert data to spatial geopandas.GeoDataFrame
    tidepoints_gdf = gpd.GeoDataFrame(
        data={"time": tidepoints_df.index, "tide_m": tidepoints_df.tide_m},
        geometry=gpd.points_from_xy(tidepoints_df.lon, tidepoints_df.lat),
        crs="EPSG:4326",
    )

    # Reproject to satellite data CRS
    tidepoints_gdf = tidepoints_gdf.to_crs(crs=ds.crs)

    # Fix time and set to index
    tidepoints_gdf["time"] = pd.to_datetime(tidepoints_gdf["time"], utc=True)
    tidepoints_gdf = tidepoints_gdf.set_index("time")

    return tidepoints_gdf


def interpolate_tide(timestep, tidepoints_gdf, method="rbf", factor=50):
    """
    Extract a subset of tide modelling point data for a given time-step,
    then interpolate these tides into the extent of the xarray dataset.

    Parameters:
    -----------
    timestep_tuple : tuple
        A tuple of x, y and time values sourced from `ds`. These values
        are used to set up the x and y grid into which tide heights for
        each timestep are interpolated. For example:
        `(ds.x.values, ds.y.values, ds.time.values)`
    tidepoints_gdf : geopandas.GeoDataFrame
        An `geopandas.GeoDataFrame` containing modelled tide heights
        with an index based on each timestep in `ds`.
    method : string, optional
        The method used to interpolate between point values. This string
        is either passed to `scipy.interpolate.griddata` (for 'linear',
        'nearest' and 'cubic' methods), or used to specify Radial Basis
        Function interpolation using `scipy.interpolate.Rbf` ('rbf').
        Defaults to 'rbf'.
    factor : int, optional
        An optional integer that can be used to subsample the spatial
        interpolation extent to obtain faster interpolation times, then
        up-sample this array back to the original dimensions of the
        data as a final step. For example, setting `factor=10` will
        interpolate ata into a grid that has one tenth of the
        resolution of `ds`. This approach will be significantly faster
        than interpolating at full resolution, but will potentially
        produce less accurate or reliable results.

    Returns:
    --------
    out_tide : xarray.DataArray
        A 2D array containing tide heights interpolated into the extent
        of the input data.
    """

    # Extract subset of observations based on timestamp of imagery
    time_string = str(timestep.time.values)[0:19].replace("T", " ")
    tidepoints_subset = tidepoints_gdf.loc[time_string]

    # Get lists of x, y and z (tide height) data to interpolate
    x_coords = tidepoints_subset.geometry.x.values.astype("float32")
    y_coords = tidepoints_subset.geometry.y.values.astype("float32")
    z_coords = tidepoints_subset.tide_m.values.astype("float32")

    # Interpolate tides into the extent of the satellite timestep
    out_tide = interpolate_2d(
        ds=timestep,
        x_coords=x_coords,
        y_coords=y_coords,
        z_coords=z_coords,
        method=method,
        factor=factor,
    )

    # Return data as a Float32 to conserve memory
    return out_tide.astype("float32")


def multiprocess_apply(ds, dim, func):
    """
    Applies a custom function along the dimension of an xarray.Dataset,
    then combines the output to match the original dataset.

    Parameters:
    -----------
    ds : xarray.Dataset
        A dataset with a dimension `dim` to apply the custom function
        along.
    dim : string
        The dimension along which the custom function will be applied.
    func : function
        The function that will be applied in parallel to each array
        along dimension `dim`. To specify custom parameters, use
        `functools.partial`.

    Returns:
    --------
    xarray.Dataset
        A concatenated dataset containing an output for each array
        along the input `dim` dimension.
    """

    pool = multiprocessing.Pool(multiprocessing.cpu_count() - 1)
    out_list = pool.map(func, iterable=[group for (i, group) in ds.groupby(dim)])
    pool.close()
    pool.join()

    # Combine to match the original dataset
    return xr.concat(out_list, dim=ds[dim])


def load_tidal_subset(year_ds, tide_cutoff_min, tide_cutoff_max):
    """
    For a given year of data, thresholds data to keep observations
    within a minimum and maximum tide height cutoff range, and load
    the data into memory.

    Parameters:
    -----------
    year_ds : xarray.Dataset
        An `xarray.Dataset` for a single epoch (typically annually)
        containing a time series of water index data (e.g. MNDWI) and
        tide heights (`tide_m`) for each pixel.
    tide_cutoff_min, tide_cutoff_max : numeric or xarray.DataArray
        Numeric values or 2D data arrays containing minimum and
        maximum tide height values used to select a subset of
        satellite observations for each individual pixel that fall
        within this range. All pixels with tide heights outside of
        this range will be set to `NaN`.

    Returns:
    --------
    year_ds : xarray.Dataset
        An in-memory `xarray.Dataset` with pixels set to `NaN` if
        they were acquired outside of the supplied tide height range.
    """

    # Determine what pixels were acquired in selected tide range, and
    # drop time-steps without any relevant pixels to reduce data to load
    tide_bool = (year_ds.tide_m >= tide_cutoff_min) & (
        year_ds.tide_m <= tide_cutoff_max
    )
    year_ds = year_ds.sel(time=tide_bool.sum(dim=["x", "y"]) > 0)

    # Apply mask, and load in corresponding tide masked data
    year_ds = year_ds.where(tide_bool)
    return year_ds.compute()


def tidal_composite(
    year_ds, label, label_dim, output_dir, output_suffix="", export_geotiff=False
):
    """
    For a given year of data, takes median, counts and standard
    deviationo of valid water index results, and optionally writes
    each water index, tide height, standard deviation and valid pixel
    counts for the time period to file as GeoTIFFs.

    Parameters:
    -----------
    year_ds : xarray.Dataset
        A tide-masked `xarray.Dataset` containing a time series of
        water index data (e.g. MNDWI) for a single epoch
        (typically annually).
    label : int, float or str
        An int, float or string used to name the output variables and
        files, For example, a year label for the data in `year_ds`.
    label_dim : str
        A string giving the name for the dimension that will be created
        to contain labels supplied to `label`.
    output_dir : str
        The directory to output files for the specific analysis.
    output_suffix : str
        An optional suffix that will be used to identify certain
        outputs. For example, outputs used for gapfilling data
        can be differentiated from other outputs by supplying
        `output_suffix='_gapfill'`. Defaults to '', which will not
        add a suffix to the output file names.
    output_geotiff : bool, optional
        Whether to export output files as GeoTIFFs. Defaults to False.

    Returns:
    --------
    median_ds : xarray.Dataset
        An in-memory `xarray.Dataset` containing output composite
        arrays with a dimension labelled by `label` and `label_dim`.
    """

    # Compute median water indices and counts of valid pixels
    median_ds = year_ds.median(dim="time", keep_attrs=True)
    median_ds["count"] = year_ds.mndwi.count(dim="time", keep_attrs=True).astype(
        "int16"
    )
    median_ds["stdev"] = year_ds.mndwi.std(dim="time", keep_attrs=True)

    # Set nodata values
    median_ds["ndwi"].attrs["nodata"] = np.nan
    median_ds["mndwi"].attrs["nodata"] = np.nan
    median_ds["tide_m"].attrs["nodata"] = np.nan
    median_ds["stdev"].attrs["nodata"] = np.nan
    median_ds["count"].attrs["nodata"] = -999

    # Write each variable to file
    if export_geotiff:
        for i in median_ds:
            write_cog(
                geo_im=median_ds[i].compute(),
                fname=f"{output_dir}/{str(label)}_{i}{output_suffix}.tif",
                overwrite=True,
            )

    # Set coordinate and dim
    median_ds = median_ds.assign_coords(**{label_dim: label}).expand_dims(label_dim)

    return median_ds


def export_annual_gapfill(ds, output_dir, tide_cutoff_min, tide_cutoff_max):
    """
    To calculate both annual median composites and three-year gapfill
    composites without having to load more than three years in memory
    at the one time, this function loops through the years in the
    dataset, progressively updating three datasets (the previous year,
    current year and subsequent year of data).

    Parameters:
    -----------
    ds : xarray.Dataset
        A tide-masked `xarray.Dataset` containing a time series of
        water index data (e.g. MNDWI).
    output_dir : str
        The directory to output files for the specific analysis.
    tide_cutoff_min, tide_cutoff_max : numeric or xarray.DataArray
        Numeric values or 2D data arrays containing minimum and
        maximum tide height values used to select a subset of
        satellite observations for each individual pixel that fall
        within this range. All pixels with tide heights outside of
        this range will be set to `NaN`.
    """

    # Create empty vars containing un-composited data from the previous,
    # current and future year. This is progressively updated to ensure that
    # no more than 3 years of data are loaded into memory at any one time
    previous_ds = None
    current_ds = None
    future_ds = None

    # Iterate through each year in the dataset, starting at one year before
    for year in np.unique(ds.time.dt.year) - 1:

        # Load data for the subsequent year
        future_ds = load_tidal_subset(
            ds.sel(time=str(year + 1)),
            tide_cutoff_min=tide_cutoff_min,
            tide_cutoff_max=tide_cutoff_max,
        )

        # If the current year var contains data, combine these observations
        # into median annual high tide composites and export GeoTIFFs
        if current_ds:

            # Generate composite
            tidal_composite(
                current_ds,
                label=year,
                label_dim="year",
                output_dir=output_dir,
                export_geotiff=True,
            )

        # If ALL of the previous, current and future year vars contain data,
        # combine these three years of observations into a single median
        # 3-year gapfill composite
        if previous_ds and current_ds and future_ds:

            # Concatenate the three years into one xarray.Dataset
            gapfill_ds = xr.concat([previous_ds, current_ds, future_ds], dim="time")

            # Generate composite
            tidal_composite(
                gapfill_ds,
                label=year,
                label_dim="year",
                output_dir=output_dir,
                output_suffix="_gapfill",
                export_geotiff=True,
            )

        # Shift all loaded data back so that we can re-use it in the next
        # iteration and not have to load the same data multiple times
        previous_ds = current_ds
        current_ds = future_ds
        future_ds = []


def generate_rasters(
    dc, config, study_area, raster_version, start_year, end_year, log=None
):
    #####################################
    # Connect to datacube, Dask cluster #
    #####################################
    if log is None:
        log = configure_logging()

    # Create local dask client for parallelisation
    client = create_local_dask_cluster(return_client=True)

    ###########################
    # Load supplementary data #
    ###########################

    # Tide points are used to model tides across the extent of the satellite data
    points_gdf = gpd.read_file(config["Input files"]["coastal_points_path"])

    # Grid cells used to process the analysis
    gridcell_gdf = (
        gpd.read_file(config["Input files"]["coastal_grid_path"])
        .to_crs(epsg=4326)
        .set_index("id")
    )
    gridcell_gdf.index = gridcell_gdf.index.astype(int).astype(str)
    gridcell_gdf = gridcell_gdf.loc[[str(study_area)]]
    log.info("Loaded tide modelling points and study area grid")

    ################
    # Loading data #
    ################

    # Create query
    geopoly = Geometry(gridcell_gdf.iloc[0].geometry, crs=gridcell_gdf.crs)
    query = {
        "geopolygon": geopoly.buffer(0.05),
        "time": (start_year, end_year),
        "dask_chunks": {"time": 1, "x": 3000, "y": 3000},
    }

    # Load virtual product
    try:
        ds = load_water_index(
            dc,
            query,
            yaml_path=config["Virtual product"]["virtual_product_path"],
            product_name=config["Virtual product"]["virtual_product_name"],
        )
    except (ValueError, IndexError):
        raise ValueError(f"No valid data found for gridcell {study_area}")
    log.info("Loaded virtual product")

    ###################
    # Tidal modelling #
    ###################

    # Model tides at each point in a provided geopandas.GeoDataFrame
    # based on all timesteps observed by Landsat. This returns a new
    # geopandas.GeoDataFrame with a "time" index (matching every time
    # step in our Landsat data), and a "tide_m" column giving the tide
    # heights at each point location at that time.
    tidepoints_gdf = model_tide_points(ds, points_gdf)

    # Test if there is data and skip rest of the analysis if not
    if tidepoints_gdf.geometry.unique().shape[0] <= 1:
        raise Exception(
            "Gridcell {study_area} has 1 or less tidal points so cannot interpolate tide data"
        )
    log.info("Modelled tide heights for each tide modelling point")

    # For each satellite timestep, spatially interpolate our modelled
    # tide height points into the spatial extent of our satellite image,
    # and add this new data as a new variable in our satellite dataset.
    # This allows each satellite pixel to be analysed and filtered
    # based on the tide height at the exact moment of each satellite
    # image acquisition.
    ds["tide_m"] = multiprocess_apply(
        ds=ds, dim="time", func=partial(interpolate_tide, tidepoints_gdf=tidepoints_gdf)
    )
    log.info("Finished spatially interpolating tide heights")

    # Based on the entire time-series of tide heights, compute the max
    # and min satellite-observed tide height for each pixel, then
    # calculate tide cutoffs used to restrict our data to satellite
    # observations centred over mid-tide (0 m Above Mean Sea Level).
    tide_cutoff_buff = (
        ds["tide_m"].max(dim="time") - ds["tide_m"].min(dim="time")
    ) * 0.25
    tide_cutoff_min = 0.0 - tide_cutoff_buff
    tide_cutoff_max = 0.0 + tide_cutoff_buff
    log.info("Calculated tide cutoffs")

    ##############################
    # Generate yearly composites #
    ##############################

    # If output folder doesn't exist, create it
    output_dir = (
        f"data/interim/raster/{raster_version}/" f"{study_area}_{raster_version}"
    )
    os.makedirs(output_dir, exist_ok=True)

    # Iterate through each year and export annual and 3-year
    # gapfill composites
    export_annual_gapfill(ds, output_dir, tide_cutoff_min, tide_cutoff_max)
    log.info("Completed writing data")

    # Close dask client
    client.close()


@click.command()
@click.option(
    "--config_path",
    type=str,
    required=True,
    help="Path to the YAML config file defining inputs to "
    "use for this analysis. These are typically located in "
    "the `dea-coastlines/configs/` directory.",
)
@click.option(
    "--study_area",
    type=str,
    required=True,
    help="A string providing a unique ID of an analysis "
    "gridcell that will be used to run the analysis. This "
    'should match a row in the "id" column of the provided '
    "analysis gridcell vector file.",
)
@click.option(
    "--raster_version",
    type=str,
    required=True,
    help="A unique string proving a name that will be used "
    "for output raster directories and files. This can be "
    "used to version different analysis outputs.",
)
@click.option(
    "--start_year",
    type=str,
    default="1987",
    help="The first year used to load data. Note that this "
    "should buffer the desired temporal extent of the "
    "analysis by one year to allow sufficient data for "
    "gapfilling low data pixels. For example, set "
    "`--start_year 1987` to extract a shoreline timeseries "
    "that commences in 1988.",
)
@click.option(
    "--end_year",
    type=str,
    default="2021",
    help="The last year used to load data. Note that this "
    "should buffer the desired temporal extent of the "
    "analysis by one year to allow sufficient data for "
    "gapfilling low data pixels. For example, set "
    "`--end_year 2021` to extract a shoreline timeseries "
    "that finishes in the year 2020.",
)
def generate_rasters_cli(config_path, study_area, raster_version, start_year, end_year):
    # Connect to datacube
    dc = datacube.Datacube(app="Coastlines")

    # Load analysis params from config file
    config = load_config(config_path=config_path)

    log = configure_logging(f"Coastlines Raster {study_area}")

    try:
        generate_rasters(
            dc,
            config,
            study_area,
            raster_version,
            start_year,
            end_year,
            log=log,
        )
    except Exception as e:
        log.exception(
            f"Failed to run process on study area {study_area} with error {e}"
        )
        sys.exit(1)


if __name__ == "__main__":
    generate_rasters_cli()
