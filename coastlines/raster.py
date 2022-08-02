#!/usr/bin/env python
# coding: utf-8

# This code conducts raster generation for DE Africa Coastlines:

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
import multiprocess
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
from dea_tools.spatial import interpolate_2d

from coastlines.utils import configure_logging, load_config

# Hide warnings
warnings.simplefilter(action="ignore", category=FutureWarning)


def hillshade(dem, elevation, azimuth, vert_exag=1, dx=30, dy=30):
    """
    Calculate hillshade from an input Digital Elevation Model
    (DEM) array and a sun elevation and azimith.

    Parameters:
    -----------
    dem : numpy.array
        A 2D Digital Elevation Model array.
    elevation : int or float
        Sun elevation (0-90, degrees up from horizontal).
    azimith : int or float
        Sun azimuth (0-360, degrees clockwise from north).
    vert_exag : int or float, optional
        The amount to exaggerate the elevation values by
        when calculating illumination. This can be used either
        to correct for differences in units between the x-y coordinate
        system and the elevation coordinate system (e.g. decimal
        degrees vs. meters) or to exaggerate or de-emphasize
        topographic effects.
    dx : int or float, optional
        The x-spacing (columns) of the input DEM. This
        is typically the spatial resolution of the DEM.
    dy : int or float, optional
        The y-spacing (rows) of the input input DEM. This
        is typically the spatial resolution of the DEM.

    Returns:
    --------
    hs : numpy.array
        A 2D hillshade array with values between 0-1, where
        0 is completely in shadow and 1 is completely
        illuminated.
    """

    from matplotlib.colors import LightSource

    hs = LightSource(azdeg=azimuth, altdeg=elevation).hillshade(
        dem, vert_exag=vert_exag, dx=dx, dy=dy
    )
    return hs


def sun_angles(dc, query):
    """
    For a given spatiotemporal query, calculate mean sun
    azimuth and elevation for each satellite observation, and
    return these as a new `xarray.Dataset` with 'sun_elevation'
    and 'sun_azimuth' variables.

    Parameters:
    -----------
    dc : datacube.Datacube object
        Datacube instance used to load data.
    query : dict
        A dictionary containing query parameters used to identify
        satellite observations and load metadata.

    Returns:
    --------
    sun_angles_ds : xarray.Dataset
        An `xarray.set` containing a 'sun_elevation' and
        'sun_azimuth' variables.
    """

    from datacube.api.query import query_group_by
    from datacube.model.utils import xr_apply

    # Identify satellite datasets and group outputs using the
    # same approach used to group satellite imagery (i.e. solar day)
    gb = query_group_by(**query)
    datasets = dc.find_datasets(**query)
    dataset_array = dc.group_datasets(datasets, gb)

    # Load and take the mean of metadata from each product
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
    terrain_shadow_ds = multiprocess_apply(
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

    return ds[["mndwi", "ndwi"]]


def model_tides(
    x,
    y,
    time,
    model="FES2014",
    directory="~/tide_models_clipped",
    epsg=4326,
    method="bilinear",
    extrapolate=True,
    cutoff=10.0,
):
    """
    Compute tides at points and times using tidal harmonics.
    If multiple x, y points are provided, tides will be
    computed for all timesteps at each point.

    This function supports any tidal model supported by
    `pyTMD`, including the FES2014 Finite Element Solution
    tide model, and the TPXO8-atlas and TPXO9-atlas-v5
    TOPEX/POSEIDON global tide models.

    This function requires access to tide model data files
    to work. These should be placed in a folder with
    subfolders matching the formats specified by `pyTMD`:
    https://pytmd.readthedocs.io/en/latest/getting_started/Getting-Started.html#directories

    For FES2014 (https://www.aviso.altimetry.fr/es/data/products/auxiliary-products/global-tide-fes/description-fes2014.html):
        - {directory}/fes2014/ocean_tide/
          {directory}/fes2014/load_tide/

    For TPXO8-atlas (https://www.tpxo.net/tpxo-products-and-registration):
        - {directory}/tpxo8_atlas/

    For TPXO9-atlas-v5 (https://www.tpxo.net/tpxo-products-and-registration):
        - {directory}/TPXO9_atlas_v5/

    This function is a minor modification of the `pyTMD`
    package's `compute_tide_corrections` function, adapted
    to process multiple timesteps for multiple input point
    locations. For more info:
    https://pytmd.readthedocs.io/en/stable/user_guide/compute_tide_corrections.html

    Parameters:
    -----------
    x, y : float or list of floats
        One or more x and y coordinates used to define
        the location at which to model tides. By default these
        coordinates should be lat/lon; use `epsg` if they
        are in a custom coordinate reference system.
    time : A datetime array or pandas.DatetimeIndex
        An array containing 'datetime64[ns]' values or a
        'pandas.DatetimeIndex' providing the times at which to
        model tides in UTC time.
    model : string
        The tide model used to model tides. Options include:
        - FES2014
        - TPXO8-atlas
        - TPXO9-atlas-v5
    directory : string
        The directory containing tide model data files. These
        data files should be stored in sub-folders for each
        model that match the structure provided by `pyTMD`:
        https://pytmd.readthedocs.io/en/latest/getting_started/Getting-Started.html#directories
        For example:
        - {directory}/fes2014/ocean_tide/
          {directory}/fes2014/load_tide/
        - {directory}/tpxo8_atlas/
        - {directory}/TPXO9_atlas_v5/
    epsg : int
        Input coordinate system for 'x' and 'y' coordinates.
        Defaults to 4326 (WGS84).
    method : string
        Method used to interpolate tidal contsituents
        from model files. Options include:
        - bilinear: quick bilinear interpolation
        - spline: scipy bivariate spline interpolation
        - linear, nearest: scipy regular grid interpolations
    extrapolate : bool
        Whether to extrapolate tides for locations outside of
        the tide modelling domain using nearest-neighbor
    cutoff : int or float
        Extrapolation cutoff in kilometers. Set to `np.inf`
        to extrapolate for all points.

    Returns
    -------
    A pandas.DataFrame containing tide heights for every
    combination of time and point coordinates.
    """
    import os
    import pyproj
    import numpy as np
    import pyTMD.time
    import pyTMD.model
    import pyTMD.utilities
    from pyTMD.calc_delta_time import calc_delta_time
    from pyTMD.infer_minor_corrections import infer_minor_corrections
    from pyTMD.predict_tide_drift import predict_tide_drift
    from pyTMD.read_tide_model import extract_tidal_constants
    from pyTMD.read_netcdf_model import extract_netcdf_constants
    from pyTMD.read_GOT_model import extract_GOT_constants
    from pyTMD.read_FES_model import extract_FES_constants

    # Check that tide directory is accessible
    try:
        os.access(directory, os.F_OK)
    except:
        raise FileNotFoundError("Invalid tide directory")

    # Get parameters for tide model
    model = pyTMD.model(directory, format="netcdf", compressed=False).elevation(model)

    # If time passed as a single Timestamp, convert to datetime64
    if isinstance(time, pd.Timestamp):
        time = time.to_datetime64()

    # Handle numeric or array inputs
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    time = np.atleast_1d(time)

    # Determine point and time counts
    assert len(x) == len(y), "x and y must be the same length"
    n_points = len(x)
    n_times = len(time)

    # Converting x,y from EPSG to latitude/longitude
    try:
        # EPSG projection code string or int
        crs1 = pyproj.CRS.from_string("epsg:{0:d}".format(int(epsg)))
    except (ValueError, pyproj.exceptions.CRSError):
        # Projection SRS string
        crs1 = pyproj.CRS.from_string(epsg)

    crs2 = pyproj.CRS.from_string("epsg:{0:d}".format(4326))
    transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
    lon, lat = transformer.transform(x.flatten(), y.flatten())

    # Assert delta time is an array and convert datetime
    time = np.atleast_1d(time)
    t = pyTMD.time.convert_datetime(time, epoch=(1992, 1, 1, 0, 0, 0)) / 86400.0

    # Delta time (TT - UT1) file
    delta_file = pyTMD.utilities.get_data_path(["data", "merged_deltat.data"])

    # Read tidal constants and interpolate to grid points
    if model.format in ("OTIS", "ATLAS"):
        amp, ph, D, c = extract_tidal_constants(
            lon,
            lat,
            model.grid_file,
            model.model_file,
            model.projection,
            TYPE=model.type,
            METHOD=method,
            EXTRAPOLATE=extrapolate,
            CUTOFF=cutoff,
            GRID=model.format,
        )
        deltat = np.zeros_like(t)

    elif model.format == "netcdf":
        amp, ph, D, c = extract_netcdf_constants(
            lon,
            lat,
            model.grid_file,
            model.model_file,
            TYPE=model.type,
            METHOD=method,
            EXTRAPOLATE=extrapolate,
            CUTOFF=cutoff,
            SCALE=model.scale,
            GZIP=model.compressed,
        )
        deltat = np.zeros_like(t)

    elif model.format == "GOT":
        amp, ph, c = extract_GOT_constants(
            lon,
            lat,
            model.model_file,
            METHOD=method,
            EXTRAPOLATE=extrapolate,
            CUTOFF=cutoff,
            SCALE=model.scale,
            GZIP=model.compressed,
        )

        # Interpolate delta times from calendar dates to tide time
        deltat = calc_delta_time(delta_file, t)

    elif model.format == "FES":
        amp, ph = extract_FES_constants(
            lon,
            lat,
            model.model_file,
            TYPE=model.type,
            VERSION=model.version,
            METHOD=method,
            EXTRAPOLATE=extrapolate,
            CUTOFF=cutoff,
            SCALE=model.scale,
            GZIP=model.compressed,
        )

        # Available model constituents
        c = model.constituents

        # Interpolate delta times from calendar dates to tide time
        deltat = calc_delta_time(delta_file, t)

    # Calculate complex phase in radians for Euler's
    cph = -1j * ph * np.pi / 180.0

    # Calculate constituent oscillation
    hc = amp * np.exp(cph)

    # Repeat constituents to length of time and number of input
    # coords before passing to `predict_tide_drift`
    t, hc, deltat = (
        np.tile(t, n_points),
        hc.repeat(n_times, axis=0),
        np.tile(deltat, n_points),
    )

    # Predict tidal elevations at time and infer minor corrections
    npts = len(t)
    tide = np.ma.zeros((npts), fill_value=np.nan)
    tide.mask = np.any(hc.mask, axis=1)
    tide.data[:] = predict_tide_drift(t, hc, c, DELTAT=deltat, CORRECTIONS=model.format)
    minor = infer_minor_corrections(t, hc, c, DELTAT=deltat, CORRECTIONS=model.format)
    tide.data[:] += minor.data[:]

    # Replace invalid values with fill value
    tide.data[tide.mask] = tide.fill_value

    # Export data as a dataframe
    return pd.DataFrame(
        {
            "time": np.tile(time, n_points),
            "x": np.repeat(x, n_times),
            "y": np.repeat(y, n_times),
            "tide_m": tide,
        }
    ).set_index("time")


def pixel_tides(
    ds,
    resample_func=None,
    times=None,
    calculate_quantiles=None,
    resolution=5000,
    buffer=12000,
    model="FES2014",
    directory="~/tide_models_clipped",
):
    """
    Obtain tide heights for each pixel in a dataset by modelling
    tides into a low-resolution grid surrounding the dataset,
    then (optionally) spatially reprojecting this low-res data back
    into the original higher resolution dataset extent and resolution.

    Parameters:
    -----------
    ds : xarray.Dataset
        A dataset whose geobox (`ds.odc.geobox`) will be used to define
        the spatial extent of the low resolution tide modelling grid.
    resample_func : function, optional
        A function to use to re-project low resolution tides back into
        `ds`'s original higher resolution grid. If you do not want low
        resolution tides to be re-projected back to higher resolution,
        set this to `None`.
    times : list or pandas.Series, optional
        By default, the function will model tides using the times
        contained in the `time` dimension of `ds`. This param can be used
        to model tides for a custom set of times instead, for example:
        `times=pd.date_range(start="2000", end="2020", freq="30min")`
    calculate_quantiles : list or np.array, optional
        Rather than returning all individual tides, low-resolution tides
        can be first aggregated using a quantile calculation by passing in
        a list or array of quantiles to compute. For example, this could
        be used to calculate the min/max tide across all times:
        `calculate_quantiles=[0.0, 1.0]`.
    resolution: int, optional
        The desired resolution of the low-resolution grid used for tide
        modelling. Defaults to 5000 for a 5,000 m grid (assuming `ds`'s
        CRS uses project/metre units).
    buffer : int, optional
        The amount by which to buffer the higher resolution grid extent
        when creating the new low resolution grid. This buffering is
        important as it ensures that ensure pixel-based tides are seamless
        across dataset boundaries. This buffer will eventually be clipped
        away when the low-resolution data is re-projected back to the
        resolution and extent of the higher resolution dataset. Defaults
        to 12,000 m to ensure that at least two 5,000 m pixels occur
        outside of the dataset bounds.
    model : string, optional
        The tide model used to model tides. Options include:
        - FES2014
        - TPXO8-atlas
        - TPXO9-atlas-v5
    directory : string, optional
        The directory containing tide model data files. These
        data files should be stored in sub-folders for each
        model that match the structure provided by `pyTMD`:
        https://pytmd.readthedocs.io/en/latest/getting_started/Getting-Started.html#directories
        For example:
        - {directory}/fes2014/ocean_tide/
          {directory}/fes2014/load_tide/
        - {directory}/tpxo8_atlas/
        - {directory}/TPXO9_atlas_v5/

    Returns:
    --------
    If `resample_func` is None:

        tides_lowres : xr.DataArray
            A low resolution data array giving either tide heights every
            timestep in `ds` (if `times` is None), tide heights at every
            time in `times` (if `times` is not None), or tide height quantiles
            for every quantile provided by `calculate_quantiles`.

    If `sample_func` is not None:

        tides_highres, tides_lowres : tuple of xr.DataArrays
            In addition to `tides_lowres` (see above), a high resolution
            array of tide heights will be generated that matches the
            exact spatial resolution and extent of `ds`. This will contain
            either tide heights every timestep in `ds` (if `times` is None),
            tide heights at every time in `times` (if `times` is not None),
            or tide height quantiles for every quantile provided by
            `calculate_quantiles`.
    """

    from odc.geo.geobox import GeoBox

    # Create a new reduced resolution (5km) tide modelling grid after
    # first buffering the grid by 12km (i.e. at least two 5km pixels)
    print("Creating reduced resolution tide modelling array")
    buffered_geobox = ds.odc.geobox.buffered(buffer)
    rescaled_geobox = GeoBox.from_bbox(
        bbox=buffered_geobox.boundingbox, resolution=resolution
    )
    rescaled_ds = odc.geo.xr.xr_zeros(rescaled_geobox)

    # Flatten grid to 1D, then add time dimension. If custom times are
    # provided use these, otherwise use times from `ds`
    time_coords = ds.coords["time"] if times is None else times
    flattened_ds = rescaled_ds.stack(z=("x", "y"))
    flattened_ds = flattened_ds.expand_dims(dim={"time": time_coords.values})

    # Model tides for each timestep and x/y grid cell
    print(f"Modelling tides using {model} tide model")
    tide_df = model_tides(
        x=flattened_ds.x,
        y=flattened_ds.y,
        time=flattened_ds.time,
        model=model,
        directory=directory,
        epsg=ds.odc.geobox.crs.epsg,
    )

    # Insert modelled tide values back into flattened array, then unstack
    # back to 3D (x, y, time)
    tides_lowres = (
        tide_df.set_index(["x", "y"], append=True)
        .to_xarray()
        .tide_m.reindex_like(rescaled_ds)
        .transpose(*list(ds.dims.keys()))
        .astype(np.float32)
    )

    # Optionally calculate and return quantiles rather than raw data
    if calculate_quantiles is not None:

        print("Computing tide quantiles")
        tides_lowres = tides_lowres.quantile(q=calculate_quantiles, dim="time")
        reproject_dim = "quantile"

    else:
        reproject_dim = "time"

    # Ensure CRS is present
    tides_lowres = tides_lowres.odc.assign_crs(ds.odc.geobox.crs)

    # Reproject each timestep into original high resolution grid
    if resample_func:

        print("Reprojecting tides into original array")
        tides_highres = parallel_apply(
            tides_lowres, reproject_dim, odc.algo.xr_reproject, ds.odc.geobox.compat
        )

        return tides_highres, tides_lowres

    else:
        print("Returning low resolution tide array")
        return tides_lowres


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
        satellite observations are approximiately centred over this
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

    from tqdm import tqdm

    with multiprocess.Pool(multiprocess.cpu_count() - 1) as pool:

        # Apply func in parallel
        to_iterate = [group for (i, group) in ds.groupby(dim)]
        out_list = list(tqdm(pool.map(func, to_iterate), total=len(to_iterate)))

    # Combine to match the original dataset
    return xr.concat(out_list, dim=ds[dim])


def parallel_apply(ds, dim, func, *args):
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
        along dimension `dim`. The first argument passed to this
        function should be the array along `dim`.
    *args :
        Any number of arguments that will be passed to `func`.
    Returns:
    --------
    xarray.Dataset
        A concatenated dataset containing an output for each array
        along the input `dim` dimension.
    """

    from concurrent.futures import ProcessPoolExecutor
    from tqdm import tqdm
    from itertools import repeat

    with ProcessPoolExecutor() as executor:

        # Apply func in parallel
        groups = [group for (i, group) in ds.groupby(dim)]
        to_iterate = (groups, *(repeat(i, len(groups)) for i in args))
        out_list = list(tqdm(executor.map(func, *to_iterate), total=len(groups)))

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
    ds["tide_m"], tides_lowres = pixel_tides(ds, resample_func=True)
    log.info(f"Study area {study_area}: Finished modelling tide heights")

    # Based on the entire time-series of tide heights, compute the max
    # and min satellite-observed tide height for each pixel, then
    # calculate tide cutoffs used to restrict our data to satellite
    # observations centred over mid-tide (0 m Above Mean Sea Level).
    tide_cutoff_min, tide_cutoff_max = tide_cutoffs(ds, tides_lowres, tide_centre=0.0)
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
    export_annual_gapfill(ds, output_dir, tide_cutoff_min, tide_cutoff_max)
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
    aws_unsigned,
    overwrite,
):
    log = configure_logging(f"Coastlines raster generation for study area {study_area}")

    # Test if study area has already been run by checking if final raster exists
    output_exists = os.path.exists(
        f"data/interim/raster/{raster_version}/{study_area}_{raster_version}/{int(end_year) - 1}_mndwi.tif"
    )

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
            log=log,
        )
    except Exception as e:
        log.exception(f"Study area {study_area}: Failed to run process with error {e}")
        sys.exit(1)


if __name__ == "__main__":
    generate_rasters_cli()
