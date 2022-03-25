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

# Standard library
import sys
import os
import datetime
import warnings
import multiprocessing
from functools import partial
from collections import Counter

# Third party
import pytz
import yaml
import click
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
from shapely.geometry import shape

# Datacube and dea-tools funcs
import datacube
import odc.algo
from datacube.utils.cog import write_cog
from datacube.utils.geometry import Geometry
from datacube.utils.masking import make_mask
from datacube.virtual import catalog_from_file, construct
from dea_tools.spatial import interpolate_2d
from dea_tools.dask import create_local_dask_cluster

# Hide warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def load_config(config_path):
    """
    Loads a YAML config file and returns data as a nested dictionary.
    """
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def load_water_index(dc, query, yaml_path, product_name='ls_nbart_mndwi'):
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
        geotransform = ds.metadata_doc['grids']['default']['transform']
        shape = ds.metadata_doc['grids']['default']['shape']
        crs = CRS(ds.metadata_doc['crs'])
        affine = Affine(geotransform[0], 0.0, geotransform[2], 0.0,
                        geotransform[4], geotransform[5])

        return GeoBox(width=shape[1],
                      height=shape[0],
                      affine=affine,
                      crs=crs)

    # Load in virtual product catalogue and select MNDWI product
    catalog = catalog_from_file(yaml_path)
    product = catalog[product_name]

    # Construct a new version of the product using most common CRS
    # Determine geobox with custom function to increase lazy loading
    # speed (will eventually be done automatically within virtual
    # products)
    with mock.patch('datacube.virtual.impl.native_geobox',
                    side_effect=custom_native_geobox):

        # Identify most common CRS
        bag = product.query(dc, **query)
        crs_list = [str(i.crs) for i in bag.contained_datasets()]
        crs_counts = Counter(crs_list)
        crs = crs_counts.most_common(1)[0][0]

        # Pass CRS to product load
        settings = dict(output_crs=crs,
                        resolution=(-30, 30),
                        align=(15, 15),
                        resampling={
                            'pixel_quality': 'nearest',
                            '*': 'cubic'
                        })
        box = product.group(bag, **settings, **query)
        ds = product.fetch(box, **settings, **query)       

    # Rechunk if smallest chunk is less than 10
    if ((len(ds.x) % 3000) <= 10) or ((len(ds.y) % 3000) <= 10):
        ds = ds.chunk({'x': 3200, 'y': 3200})

    # Identify pixels that are either cloud, cloud shadow or nodata
    nodata = make_mask(ds['pixel_quality'], nodata=True)
    mask = (make_mask(ds['pixel_quality'],
                      cloud='high_confidence') |
            make_mask(ds['pixel_quality'],
                      cloud_shadow='high_confidence') | nodata)

    # Apply opening to remove long narrow false positive clouds along
    # the coastline, then dilate to restore cloud edges
    mask_cleaned = odc.algo.mask_cleanup(mask,
                                         mask_filters=[('opening', 20), 
                                                       ('dilation', 5)])
    ds = ds.where(~mask_cleaned & ~nodata)
    
    # Mask any invalid pixel values outside of 0 and 1
    green_bool = (ds.green >= 0) & (ds.green <= 1)
    swir_bool = (ds.swir_1 >= 0) & (ds.swir_1 <= 1)
    ds = ds.where(green_bool & swir_bool)
    
    # Compute MNDWI
    ds[['mndwi']] = (ds.green - ds.swir_1) / (ds.green + ds.swir_1)

    return ds[['mndwi']]


def model_tides(x,
                y,
                time,
                model="FES2014",
                directory="tide_models",
                epsg=4326,
                method="bilinear",
                extrapolate=True,
                cutoff=10.0):
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
    time : A datetime array
        An array containing 'datetime64[ns]' values providing
        the times at which to model tides in UTC time.
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

    # Determine point and time input counts
    n_points = len(x)
    n_times = len(time)

    # Verify coordinate dimension shapes
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)

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
    return pd.DataFrame({"time": np.tile(time, n_points),
                         "x": np.repeat(x, n_times),
                         "y": np.repeat(y, n_times),
                         "tide_m": tide}).set_index("time")


def model_tide_points(ds, 
                      points_gdf, 
                      extent_buffer=0.05,
                      tide_model="",
                      directory="../tide_models",
                      
                     ):
    """
    Takes an xarray.Dataset (`ds`), extracts a subset of tide modelling 
    points from a geopandas.GeoDataFrame based on`ds`'s extent, then 
    uses the `model_tides` function based on `pyTMD` to model tide 
    heights for every point at every time step in `ds`.
    
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
        model tides.
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
    observed_datetimes = ds.time.data.astype('M8[s]').astype('O').tolist()

    # Model tides for each coordinate and time
    tidepoints_df = model_tides(
        x=x_vals,
        y=y_vals,
        time=ds.time.values,
        directory=directory,
        model='FES2014')
    
    # Convert data to spatial geopandas.GeoDataFrame
    tidepoints_gdf = gpd.GeoDataFrame(data={'time': tidepoints_df.index,
                                            'tide_m': tidepoints_df.tide_m},
                                      geometry=gpd.points_from_xy(
                                          tidepoints_df.x,
                                          tidepoints_df.y),
                                      crs='EPSG:4326')

    # Reproject to satellite data CRS
    tidepoints_gdf = tidepoints_gdf.to_crs(crs=ds.crs)

    # Fix time and set to index
    tidepoints_gdf['time'] = pd.to_datetime(tidepoints_gdf['time'], utc=True)
    tidepoints_gdf = tidepoints_gdf.set_index('time')

    return tidepoints_gdf


def interpolate_tide(timestep, tidepoints_gdf, method='rbf', factor=50):
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
    time_string = str(timestep.time.values)[0:19].replace('T', ' ')
    tidepoints_subset = tidepoints_gdf.loc[time_string]
    print(f'{time_string:<80}', end='\r')

    # Get lists of x, y and z (tide height) data to interpolate
    x_coords = tidepoints_subset.geometry.x.values.astype('float32')
    y_coords = tidepoints_subset.geometry.y.values.astype('float32')
    z_coords = tidepoints_subset.tide_m.values.astype('float32')

    # Interpolate tides into the extent of the satellite timestep
    out_tide = interpolate_2d(ds=timestep,
                              x_coords=x_coords,
                              y_coords=y_coords,
                              z_coords=z_coords,
                              method=method,
                              factor=factor)

    # Return data as a Float32 to conserve memory
    return out_tide.astype('float32')


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
    print(f'Parallelising {multiprocessing.cpu_count() - 1} processes')
    out_list = pool.map(func,
                        iterable=[group for (i, group) in ds.groupby(dim)])
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

    # Print status
    year = year_ds.time[0].dt.year.item()
    print(f'Processing {year:<80}', end='\r')

    # Determine what pixels were acquired in selected tide range, and
    # drop time-steps without any relevant pixels to reduce data to load
    tide_bool = ((year_ds.tide_m >= tide_cutoff_min) &
                 (year_ds.tide_m <= tide_cutoff_max))
    year_ds = year_ds.sel(time=tide_bool.sum(dim=['x', 'y']) > 0)

    # Apply mask, and load in corresponding tide masked data
    year_ds = year_ds.where(tide_bool)
    return year_ds.compute()


def tidal_composite(year_ds,
                    label,
                    label_dim,
                    output_dir,
                    output_suffix='',
                    export_geotiff=False):
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
    median_ds = year_ds.median(dim='time', keep_attrs=True)
    median_ds['count'] = (year_ds.mndwi.count(dim='time',
                                              keep_attrs=True).astype('int16'))
    median_ds['stdev'] = year_ds.mndwi.std(dim='time', keep_attrs=True)

    # Set nodata values
    median_ds['mndwi'].attrs['nodata'] = np.nan
    median_ds['tide_m'].attrs['nodata'] = np.nan
    median_ds['stdev'].attrs['nodata'] = np.nan
    median_ds['count'].attrs['nodata'] = -999

    # Write each variable to file
    if export_geotiff:
        for i in median_ds:
            write_cog(geo_im=median_ds[i].compute(),
                      fname=f'{output_dir}/{str(label)}_{i}{output_suffix}.tif',
                      overwrite=True)

    # Set coordinate and dim
    median_ds = (median_ds.assign_coords(**{
        label_dim: label
    }).expand_dims(label_dim))

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
        future_ds = load_tidal_subset(ds.sel(time=str(year + 1)),
                                      tide_cutoff_min=tide_cutoff_min,
                                      tide_cutoff_max=tide_cutoff_max)

        # If the current year var contains data, combine these observations
        # into median annual high tide composites and export GeoTIFFs
        if current_ds:

            # Generate composite
            tidal_composite(current_ds,
                            label=year,
                            label_dim='year',
                            output_dir=output_dir,
                            export_geotiff=True)

        # If ALL of the previous, current and future year vars contain data,
        # combine these three years of observations into a single median
        # 3-year gapfill composite
        if previous_ds and current_ds and future_ds:

            # Concatenate the three years into one xarray.Dataset
            gapfill_ds = xr.concat([previous_ds, current_ds, future_ds],
                                   dim='time')

            # Generate composite
            tidal_composite(gapfill_ds,
                            label=year,
                            label_dim='year',
                            output_dir=output_dir,
                            output_suffix='_gapfill',
                            export_geotiff=True)

        # Shift all loaded data back so that we can re-use it in the next
        # iteration and not have to load the same data multiple times
        previous_ds = current_ds
        current_ds = future_ds
        future_ds = []


@click.command()
@click.option('--config_path',
              type=str,
              required=True,
              help='Path to the YAML config file defining inputs to '
              'use for this analysis. These are typically located in '
              'the `dea-coastlines/configs/` directory.')
@click.option('--study_area',
              type=str,
              required=True,
              help='A string providing a unique ID of an analysis '
              'gridcell that will be used to run the analysis. This '
              'should match a row in the "id" column of the provided '
              'analysis gridcell vector file.')
@click.option('--raster_version',
              type=str,
              required=True,
              help='A unique string proving a name that will be used '
              'for output raster directories and files. This can be '
              'used to version different analysis outputs.')
@click.option('--start_year',
              type=str,
              default='1987',
              help='The first year used to load data. Note that this '
              'should buffer the desired temporal extent of the '
              'analysis by one year to allow sufficient data for '
              'gapfilling low data pixels. For example, set '
              '`--start_year 1987` to extract a shoreline timeseries '
              'that commences in 1988.')
@click.option('--end_year',
              type=str,
              default='2021',
              help='The last year used to load data. Note that this '
              'should buffer the desired temporal extent of the '
              'analysis by one year to allow sufficient data for '
              'gapfilling low data pixels. For example, set '
              '`--end_year 2021` to extract a shoreline timeseries '
              'that finishes in the year 2020.')
def generate_rasters(config_path, study_area, raster_version, start_year,
                     end_year):

    #####################################
    # Connect to datacube, Dask cluster #
    #####################################

    # Connect to datacube
    dc = datacube.Datacube(app='DEACoastlines')

    # Create local dask client for parallelisation
    client = create_local_dask_cluster(return_client=True)

    # Load analysis params from config file
    config = load_config(config_path=config_path)

    ###########################
    # Load supplementary data #
    ###########################

    # Tide points are used to model tides across the extent of the satellite data
    points_gdf = gpd.read_file(config['Input files']['coastal_points_path'])

    # Albers grid cells used to process the analysis
    gridcell_gdf = (gpd.read_file(
        config['Input files']['coastal_grid_path']).to_crs(
            epsg=4326).set_index('id'))
    gridcell_gdf.index = gridcell_gdf.index.astype(int).astype(str)
    gridcell_gdf = gridcell_gdf.loc[[str(study_area)]]

    ################
    # Loading data #
    ################

    # Create query
    geopoly = Geometry(gridcell_gdf.iloc[0].geometry, crs=gridcell_gdf.crs)
    query = {
        'geopolygon': geopoly.buffer(0.05),
        'time': (start_year, end_year),
        'dask_chunks': {
            'time': 1,
            'x': 3000,
            'y': 3000
        }
    }

    # Load virtual product
    try:
        ds = load_water_index(
            dc,
            query,
            yaml_path=config['Virtual product']['virtual_product_path'],
            product_name=config['Virtual product']['virtual_product_name'])
    except ValueError:
            print(f'WARNING: No valid data found for gridcell {study_area}')
            sys.exit(0)     

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
        print(f'WARNING: Gridcell {study_area} has 1 or less tidal points so cannot interpolate tide data')
        sys.exit(0)

    # For each satellite timestep, spatially interpolate our modelled
    # tide height points into the spatial extent of our satellite image,
    # and add this new data as a new variable in our satellite dataset.
    # This allows each satellite pixel to be analysed and filtered
    # based on the tide height at the exact moment of each satellite
    # image acquisition.
    ds['tide_m'] = multiprocess_apply(ds=ds,
                                      dim='time',
                                      func=partial(
                                          interpolate_tide,
                                          tidepoints_gdf=tidepoints_gdf))

    # Based on the entire time-series of tide heights, compute the max
    # and min satellite-observed tide height for each pixel, then
    # calculate tide cutoffs used to restrict our data to satellite
    # observations centred over mid-tide (0 m Above Mean Sea Level).
    tide_cutoff_buff = (
        (ds['tide_m'].max(dim='time') - ds['tide_m'].min(dim='time')) * 0.25)
    tide_cutoff_min = 0.0 - tide_cutoff_buff
    tide_cutoff_max = 0.0 + tide_cutoff_buff

    ##############################
    # Generate yearly composites #
    ##############################

    # If output folder doesn't exist, create it
    output_dir = f'data/interim/raster/{raster_version}/' \
                 f'{study_area}_{raster_version}'
    os.makedirs(output_dir, exist_ok=True)

    # Iterate through each year and export annual and 3-year
    # gapfill composites
    export_annual_gapfill(ds, output_dir, tide_cutoff_min, tide_cutoff_max)

    # Close dask client
    client.close()


if __name__ == "__main__":
    generate_rasters()
