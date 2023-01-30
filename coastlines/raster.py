#!/usr/bin/env python
# coding: utf-8

# This code conducts raster generation for DEA Coastlines:

#     * Load stack of all available Landsat 5, 7, 8 and 9 satellite imagery
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

import os
import sys
import warnings
from functools import partial
from collections import Counter

import pytz
import dask
import click
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
from affine import Affine
from shapely.geometry import shape

import datacube
import odc.algo
import odc.geo.xr
from datacube.utils.aws import configure_s3_access
from datacube.utils.cog import write_cog
from datacube.utils.geometry import CRS, GeoBox, Geometry
from datacube.utils.masking import make_mask
from datacube.virtual import catalog_from_file

from dea_tools.dask import create_local_dask_cluster
from dea_tools.spatial import hillshade, sun_angles
from dea_tools.coastal import model_tides, pixel_tides
from dea_tools.datahandling import parallel_apply

from coastlines.utils import configure_logging, load_config

# Hide warnings
warnings.simplefilter(action="ignore", category=FutureWarning)


def terrain_shadow(ds, dem, threshold=0.5, radius=1):
    """
    Calculates a custom terrain shadow mask that can be used
    to remove shadow-affected satellite pixels.

    This is achieved by by first calculating hillshade using
    'sun_elevation' and 'sun_azimuth' metadata, thresholding this
    hillshade to identify shadows, then cleaning and dilating
    to return clean areas of terrain shadow.

    Parameters:
    -----------
    ds : xarray.Dataset
        An `xarray.Dataset` containing data variables named
        'sun_elevation' and 'sun_azimuth'.
    dem : numpy.array
        A 2D Digital Elevation Model array.
    threshold : float, optional
        An illumination value below which hillshaded cells should
        be considered terrain shadow. Defaults to 0.5.
    radius : int, optional
        The disk radius to use when applying morphological opening
        (to clean up small noisy shadows) and dilation (to
        mask out shadow edges. Defaults to a 1.

    Returns:
    --------
    xarray.DataArray
        An `xarray.DataArray` containing boolean True values where
        a pixel contains terrain shadow, and False if not.
    """

    from skimage.morphology import binary_dilation, binary_opening, disk

    hs = hillshade(dem, ds.sun_elevation, ds.sun_azimuth)
    hs = hs < threshold
    hs = binary_opening(hs, disk(radius))
    hs = binary_dilation(hs, disk(radius))

    return xr.DataArray(hs, dims=["y", "x"])


def terrain_shadow_masking(dc, query, ds, dem_product="dem_cop_30"):
    """
    Use a Digital Elevation Model to calculate and apply a terrain
    shadow mask to set all satellite pixels to `NaN` if they are
    affected by terrain shadow. This helps to remove noisy shorelines
    along coastal cliffs and steep coastal terrain.

    Parameters:
    -----------
    dc : datacube.Datacube object
        Datacube instance used to load DEM data and sun angle metadata.
    query : dict
        A dictionary containing query parameters used to identify
        satellite observations and load metadata.
    ds : xarray.Dataset
        An `xarray.Dataset` containing a time series of water index
        data (e.g. MNDWI) that will be terrain shadow masked.
    dem_product : string, optional
        A string giving the name of the DEM product to use for terrain
        masking. This must match the name of a product in the datacube.
        The DEM should contain a variable named 'elevation'.
        Defaults to "dem_cop_30".

    Returns:
    --------
    xarray.Dataset
        An `xarray.Dataset` containing a time series of water index
        data (e.g. MNDWI) with terrain shadowed pixels set to `NaN`.
    """

    # Compute solar angles for all satellite image timesteps
    query_subset = {k: v for k, v in query.items() if k not in ["dask_chunks"]}
    query_subset.update(
        product=["ls5_sr", "ls7_sr", "ls8_sr", "ls9_sr"],
        collection_category="T1",
        group_by="solar_day",
    )
    sun_angles_ds = sun_angles(dc, query_subset)

    # Load DEM into satellite data geobox
    dem_ds = dc.load(product=dem_product, like=ds.geobox, resampling="cubic").squeeze(
        "time", drop=True
    )
    dem_ds = dem_ds.where(dem_ds.elevation >= 0)

    # Identify terrain shadow across all timesteps
    terrain_shadow_ds = parallel_apply(
        sun_angles_ds,
        dim="time",
        func=partial(terrain_shadow, dem=dem_ds.elevation.values),
    )

    # Remove terrain shadow pixels from satellite data
    return ds.where(~terrain_shadow_ds)


def load_water_index(
    dc, query, yaml_path, product_name="ls_nbart_ndwi", mask_terrain_shadow=True
):
    """
    This function uses virtual products to load Landsat 5, 7, 8 and 9 data,
    calculate a custom remote sensing index, and return the data as a
    single xarray.Dataset.

    To minimise resampling effects and maintain the highest data
    fidelity required for subpixel shoreline extraction, this workflow
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
    mask_terrain_shadow : bool, optional
        Whether to use hillshading to mask out pixels potentially
        affected by terrain shadow. This can significantly improve
        shoreline mapping in areas of coastal cliffs or steep coastal
        topography. Defaults to True.

    Returns:
    --------
    ds : xarray.Dataset
        An `xarray.Dataset` containing a time series of water index
        data (e.g. MNDWI) for the provided datacube query
    """

    # Load in virtual product catalogue and select MNDWI product
    catalog = catalog_from_file(yaml_path)
    product = catalog[product_name]

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
        skip_broken_datasets=True,
        resampling={
            "oa_nbart_contiguity": "nearest",
            "oa_fmask": "nearest",
            "*": "cubic",
        },
    )
    box = product.group(bag, **settings, **query)
    ds = product.fetch(box, **settings, **query)

    # Rechunk if smallest chunk is less than 10
    if ((len(ds.x) % 2048) <= 10) or ((len(ds.y) % 2048) <= 10):
        ds = ds.chunk({"x": 3000, "y": 3000})

    # Extract boolean mask
    mask = odc.algo.enum_to_bool(
        ds.cloud_mask, categories=["nodata", "cloud", "shadow", "snow"]
    )

    # Close mask to remove small holes in cloud, open mask to
    # remove narrow false positive cloud, then dilate
    mask_cleaned = odc.algo.mask_cleanup(
        mask=mask, mask_filters=[("closing", 2), ("opening", 10), ("dilation", 5)]
    )

    # Add new mask to nodata pixels
    ds = odc.algo.erase_bad(ds, mask_cleaned, nodata=np.nan)

    # Apply terrain mask to remove deep shadows that can be
    # be mistaken for water
    if mask_terrain_shadow:
        ds = terrain_shadow_masking(dc, query, ds, dem_product="dem_cop_30")

    return ds[["mndwi"]]


def tide_cutoffs(ds, tides_lowres, tide_centre=0.0, resampling="bilinear"):
    """
    Based on the entire time-series of tide heights, compute the max
    and min satellite-observed tide height for each pixel, then
    calculate tide cutoffs used to restrict our data to satellite
    observations centred over mid-tide (0 m Above Mean Sea Level).

    These tide cutoffs are spatially interpolated into the extent of
    the input satellite imagery so they can be used to mask out low
    and high tide satellite pixels.

    Parameters:
    -----------
    ds : xarray.Dataset
        An `xarray.Dataset` containing a time series of water index
        data (e.g. MNDWI) for the provided datacube query. This is
        used to define the spatial extents into which tide height
        cutoffs will be interpolated.
    tides_lowres : xarray.Dataset
        A low-res `xarray.Dataset` containing tide heights for each
        timestep in `ds`, as produced by the `pixel_tides` function.
    tide_centre : float, optional
        The central tide height used to compute the min and max
        tide height cutoffs. Tide heights will be masked so all
        satellite observations are approximately centred over this
        value. The default is 0.0 which represents 0 m Above Mean
        Sea Level.
    resampling : string, optional
        The resampling method used when reprojecting low resolution
        tides to higher resolution. Defaults to "bilinear".

    Returns:
    --------
    tide_cutoff_min, tide_cutoff_max : xarray.DataArray
        2D arrays containing tide height cutoff values interpolated
        into the extent of `ds`.
    """

    # Calculate min and max tides
    tide_min = tides_lowres.min(dim="time")
    tide_max = tides_lowres.max(dim="time")

    # Identify cutoffs
    tide_cutoff_buffer = (tide_max - tide_min) * 0.25
    tide_cutoff_min = tide_centre - tide_cutoff_buffer
    tide_cutoff_max = tide_centre + tide_cutoff_buffer

    # Reproject into original geobox
    tide_cutoff_min = tide_cutoff_min.odc.reproject(
        ds.odc.geobox, resampling=resampling
    )
    tide_cutoff_max = tide_cutoff_max.odc.reproject(
        ds.odc.geobox, resampling=resampling
    )

    return tide_cutoff_min, tide_cutoff_max


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
    deviation of valid water index results, and optionally writes
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
    median_ds["stdev"] = year_ds.mndwi.std(dim="time", keep_attrs=True)
    median_ds["count"] = year_ds.mndwi.count(dim="time", keep_attrs=True).astype(
        "int16"
    )

    # Set nodata values, using np.nan for floats and -999 for ints
    for var_name, var in median_ds.data_vars.items():
        median_ds[var_name].attrs["nodata"] = -999 if var.dtype == "int16" else np.nan

    # Load data into memory
    median_ds.load()

    # Write each variable to file
    if export_geotiff:
        for i in median_ds:
            write_cog(
                geo_im=median_ds[i],
                fname=f"{output_dir}/{str(label)}_{i}{output_suffix}.tif",
                overwrite=True,
            )

    # Set coordinate and dim
    median_ds = median_ds.assign_coords(**{label_dim: label}).expand_dims(label_dim)

    return median_ds


def export_annual_gapfill(
    ds, output_dir, tide_cutoff_min, tide_cutoff_max, start_year, end_year
):
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
    start_year, end year : int
        The first and last years you wish to export annual median
        composites and three-year gapfill composites for.
    """

    # Create empty vars containing un-composited data from the previous,
    # current and future year. This is progressively updated to ensure that
    # no more than 3 years of data are loaded into memory at any one time
    previous_ds = None
    current_ds = None
    future_ds = None

    # Iterate through each year in the dataset, starting at one year before
    for year in np.arange(start_year - 2, end_year + 1):

        try:

            # Load data for the subsequent year; drop tide variable as
            # we do not need to create annual composites from this data
            future_ds = load_tidal_subset(
                ds.sel(time=str(year + 1)),
                tide_cutoff_min=tide_cutoff_min,
                tide_cutoff_max=tide_cutoff_max,
            ).drop_vars("tide_m")

        except KeyError:

            # Create empty array if error is raised due to no data being
            # available for time period
            future_ds = xr.DataArray(
                data=np.empty(shape=(0, len(ds.y), len(ds.x))),
                dims=["time", "y", "x"],
                coords={"x": ds.x, "y": ds.y, "time": []},
                name="mndwi",
            ).to_dataset()

        # If the current year var contains data, combine these observations
        # into annual median composites and export GeoTIFFs
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
    dc, config, study_area, raster_version, start_year, end_year, tide_centre, log=None
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

    # Grid cells used to process the analysis
    gridcell_gdf = (
        gpd.read_file(config["Input files"]["grid_path"])
        .to_crs(epsg=4326)
        .set_index("id")
    )
    gridcell_gdf.index = gridcell_gdf.index.astype(int).astype(str)
    gridcell_gdf = gridcell_gdf.loc[[str(study_area)]]
    log.info(f"Study area {study_area}: Loaded study area grid")

    ################
    # Loading data #
    ################

    # Create query; start year and end year are buffered by one year
    # on either side to facilitate gapfilling low data observations
    geopoly = Geometry(gridcell_gdf.iloc[0].geometry, crs=gridcell_gdf.crs)
    query = {
        "geopolygon": geopoly.buffer(0.05),
        "time": (str(start_year - 1), str(end_year + 1)),
        "dask_chunks": {"time": 1, "x": 2048, "y": 2048},
    }

    # Load virtual product
    try:
        ds = load_water_index(
            dc,
            query,
            yaml_path=config["Virtual product"]["virtual_product_path"],
            product_name=config["Virtual product"]["virtual_product_name"],
            mask_terrain_shadow=False,
        )
    except (ValueError, IndexError):
        raise ValueError(f"Study area {study_area}: No valid data found")
    log.info(f"Study area {study_area}: Loaded virtual product")

    ###################
    # Tidal modelling #
    ###################

    # For each satellite timestep, model tide heights into a low-resolution
    # 5 x 5 km grid (matching resolution of the FES2014 tidal model), then
    # reproject modelled tides into the spatial extent of our satellite image.
    # Add  this new data as a new variable in our satellite dataset to allow
    # each satellite pixel to be analysed and filtered/masked based on the
    # tide height at the exact moment of satellite image acquisition.
    try: 
        ds["tide_m"], tides_lowres = pixel_tides(ds, resample=True, directory='blah')
        log.info(f"Study area {study_area}: Finished modelling tide heights")
        
    except FileNotFoundError:
    
        log.exception(f"Study area {study_area}: Unable to access tide modelling files")
        sys.exit(2)

    # Based on the entire time-series of tide heights, compute the max
    # and min satellite-observed tide height for each pixel, then
    # calculate tide cutoffs used to restrict our data to satellite
    # observations centred over mid-tide (0 m Above Mean Sea Level).
    tide_cutoff_min, tide_cutoff_max = tide_cutoffs(
        ds, tides_lowres, tide_centre=tide_centre
    )
    log.info(
        f"Study area {study_area}: Calculating low and high tide cutoffs for each pixel"
    )

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
    log.info(f"Study area {study_area}: Started exporting raster data")
    export_annual_gapfill(
        ds, output_dir, tide_cutoff_min, tide_cutoff_max, start_year, end_year
    )
    log.info(f"Study area {study_area}: Completed exporting raster data")

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
    type=int,
    default=2000,
    help="The first annual shoreline you wish to be included "
    "in the final outputs. To allow low data pixels to be "
    "gapfilled with additional satellite data from neighbouring "
    "years, the full timeseries of satellite data loaded in this "
    "step will include one additional year of preceding satellite data "
    "(i.e. if `--start_year 2000`, satellite data from 1999 onward "
    "will be loaded for gapfilling purposes). Because of this, we "
    "recommend that at least one year of satellite data exists in "
    "your datacube prior to `--start_year`.",
)
@click.option(
    "--end_year",
    type=int,
    default=2020,
    help="The final annual shoreline you wish to be included "
    "in the final outputs. To allow low data pixels to be "
    "gapfilled with additional satellite data from neighbouring "
    "years, the full timeseries of satellite data loaded in this "
    "step will include one additional year of ensuing satellite data "
    "(i.e. if `--end_year 2020`, satellite data up to and including "
    "2021 will be loaded for gapfilling purposes). Because of this, we "
    "recommend that at least one year of satellite data exists in your "
    "datacube after `--end_year`.",
)
@click.option(
    "--tide_centre",
    type=float,
    default=0.0,
    help="The central tide height used to compute the min and max tide "
    "height cutoffs. Tide heights will be masked so all satellite "
    "observations are approximately centred over this value. The "
    "default is 0.0 which represents 0 m Above Mean Sea Level.",
)
@click.option(
    "--aws_unsigned/--no-aws_unsigned",
    type=bool,
    default=True,
    help="Whether to use sign AWS requests for S3 access",
)
@click.option(
    "--overwrite/--no-overwrite",
    type=bool,
    default=True,
    help="Whether to overwrite tiles with existing outputs, "
    "or skip these tiles entirely.",
)
def generate_rasters_cli(
    config_path,
    study_area,
    raster_version,
    start_year,
    end_year,
    tide_centre,
    aws_unsigned,
    overwrite,
):

    log = configure_logging(f"Coastlines raster generation for study area {study_area}")

    # Test if study area has already been run by checking if run status file exists
    run_status_file = f"data/interim/raster/{raster_version}/{study_area}_{raster_version}/run_completed"
    output_exists = os.path.exists(run_status_file)

    # Skip if outputs exist but overwrite is False
    if output_exists and not overwrite:
        log.info(
            f"Study area {study_area}: Data exists but overwrite set to False; skipping."
        )
        sys.exit(0)

    # Connect to datacube
    dc = datacube.Datacube(app="Coastlines")

    # Load analysis params from config file
    config = load_config(config_path=config_path)

    # Do an opinionated configuration of S3
    configure_s3_access(cloud_defaults=True, aws_unsigned=aws_unsigned)

    try:
        generate_rasters(
            dc,
            config,
            study_area,
            raster_version,
            start_year,
            end_year,
            tide_centre,
            log=log,
        )

        # Create blank run status file to indicate run completion
        with open(
            run_status_file,
            mode="w",
        ):
            pass

    except Exception as e:
        log.exception(f"Study area {study_area}: Failed to run process with error {e}")
        sys.exit(1)


if __name__ == "__main__":
    generate_rasters_cli()
