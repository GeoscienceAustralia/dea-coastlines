#!/usr/bin/env python
# coding: utf-8

# This code conducts raster generation for DEA Coastlines:

#     * Load stack of all available Landsat 5, 7 and 8 satellite imagery
#       for a location using ODC Virtual Products
#     * Convert each satellite image into a remote sensing water index 
#       (MNDWI)
#     * For each satellite image, model ocean tides into a 2 x 2 km grid
#       based on exact time of image acquisition
#     * Interpolate tide heights into spatial extent of image stack
#     * Mask out high and low tide pixels by removing all observations 
#       acquired outside of 50 percent of the observed tidal range 
#       centered over mean sea level
#     * Combine tidally-masked data into annual median composites from 
#       1988 to the present representing the coastline at approximately 
#       mean sea level (0 m AHD)
#
# Compatability:
#
#     module use /g/data/v10/public/modules/modulefiles
#     module load dea/20200713
#     pip install --user ruptures
#     pip install --user git+https://github.com/mattijn/topojson/
#     pip install --user --upgrade --extra-index-url="https://packages.dea.ga.gov.au" odc-algo
#     pip install --upgrade dask==2021.1.1 


import os
import sys
import mock
import otps
import datacube
import datetime
import odc.algo
import multiprocessing
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import scipy.interpolate
from affine import Affine
from functools import partial
from collections import Counter
from shapely.geometry import shape
from datacube.utils.cog import write_cog
from datacube.utils.dask import start_local_dask
from datacube.utils.geometry import GeoBox, Geometry, CRS, gbox
from datacube.virtual import catalog_from_file, construct
from dea_tools.spatial import interpolate_2d

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

start_time = datetime.datetime.now()


def load_mndwi(dc, 
               query,
               yaml_path, 
               product_name='ls_nbart_mndwi'):
    """
    This function uses virtual products to load data from GA Collection 
    3 Landsat 5, 7 and 8, calculate custom remote sensing indices, and 
    return the data as a single xarray.Dataset.
    
    To minimise resampling effects and maintain the highest data 
    fidelity required for subpixel coastline extraction, this workflow 
    applies masking and index calculation at native resolution, and 
    only re-projects to the most common CRS for the query using average 
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
    
    def custom_native_geobox(ds, measurements=None, basis=None):
        """
        Obtains native geobox info from dataset metadata
        """
        geotransform = ds.metadata_doc['grids']['default']['transform']
        shape = ds.metadata_doc['grids']['default']['shape']
        crs = CRS(ds.metadata_doc['crs'])
        affine = Affine(geotransform[0], 0.0, 
                        geotransform[2], 0.0, 
                        geotransform[4], geotransform[5])
        
        return GeoBox(width=shape[1], height=shape[0], 
                        affine=affine, crs=crs)  

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
                        resampling={'fmask': 'nearest', 
                                    'oa_fmask': 'nearest', 
                                    'nbart_contiguity': 'nearest',
                                    'oa_nbart_contiguity': 'nearest',
                                    '*': 'cubic'})
        box = product.group(bag, **settings, **query)
        ds = product.fetch(box, **settings, **query)   
        
    # Rechunk if smallest chunk is less than 10
    if ((len(ds.x) % 2000) <= 10) or ((len(ds.y) % 2000) <= 10):
        ds = ds.chunk({'x': 3200, 'y': 3200})

    # Extract boolean mask
    mask = odc.algo.enum_to_bool(ds.fmask,
                                 categories=['nodata', 'cloud', 'shadow', 'snow'])

    # Close mask to remove small holes in cloud, open mask to
    # remove narrow false positive cloud, then dilate
    mask_cleaned = odc.algo.mask_cleanup(mask=mask,
                                         mask_filters=[('closing', 2), 
                                                       ('opening', 10), 
                                                       ('dilation', 5)])

    # Add new mask as nodata pixels
    ds = odc.algo.erase_bad(ds, mask_cleaned, nodata=np.nan)
        
    return ds.drop('fmask')


def model_tides(ds, points_gdf, extent_buffer=0.05):
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
    
    # Obtain extent of loaded data, and f to ensure that tides are
    # modelled reliably and comparably across grid tiles
    ds_extent = shape(ds.geobox.geographic_extent.json)
    buffered = ds_extent.buffer(extent_buffer)
    subset_gdf = points_gdf[points_gdf.geometry.intersects(buffered)]

    # Extract lon, lat from tides, and time from satellite data
    x_vals = subset_gdf.geometry.centroid.x
    y_vals = subset_gdf.geometry.centroid.y
    observed_datetimes = ds.time.data.astype('M8[s]').astype('O').tolist()

    # Create list of lat/lon/time scenarios to model
    observed_timepoints = [otps.TimePoint(lon, lat, date) 
                           for date in observed_datetimes
                           for lon, lat in zip(x_vals, y_vals)]

    # Model tides for each scenario
    observed_predictedtides = otps.predict_tide(observed_timepoints)

    # Output results into pandas.DataFrame
    tidepoints_df = pd.DataFrame([(i.timepoint.timestamp, 
                                   i.timepoint.lon, 
                                   i.timepoint.lat, 
                                   i.tide_m) for i in observed_predictedtides], 
                                 columns=['time', 'lon', 'lat', 'tide_m']) 

    # Convert data to spatial geopandas.GeoDataFrame
    tidepoints_gdf = gpd.GeoDataFrame(data={'time': tidepoints_df.time, 
                                            'tide_m': tidepoints_df.tide_m}, 
                                      geometry=gpd.points_from_xy(tidepoints_df.lon, 
                                                                  tidepoints_df.lat), 
                                      crs={'init': 'EPSG:4326'})

    # Reproject to satellite data CRS
    tidepoints_gdf = tidepoints_gdf.to_crs(crs={'init': ds.crs})

    # Fix time and set to index
    tidepoints_gdf['time'] = pd.to_datetime(tidepoints_gdf['time'], utc=True)
    tidepoints_gdf = tidepoints_gdf.set_index('time')
    
    return tidepoints_gdf


def interpolate_tide(timestep, 
                     tidepoints_gdf, 
                     method='rbf', 
                     factor=20):    
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
        along dimension `dim`.
    
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
    print(f'Processing {year}')
    
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
    median_ds['count'] = (year_ds.mndwi
                          .count(dim='time', keep_attrs=True)
                          .astype('int16'))
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
    median_ds = (median_ds
                 .assign_coords(**{label_dim: label})
                 .expand_dims(label_dim)) 
        
    return median_ds


def export_annual_gapfill(ds, 
                          output_dir, 
                          tide_cutoff_min, 
                          tide_cutoff_max):
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

        

def main(argv=None):
    
    #########
    # Setup #
    #########

    if argv is None:

        argv = sys.argv
        print(sys.argv)

    # If no user arguments provided
    if len(argv) < 3:

        str_usage = "You must specify a study area ID and raster_version"
        print(str_usage)
        sys.exit()
        
    # Set study area and name for analysis
    study_area = int(argv[1])
    raster_version = str(argv[2])
    
    # Use vector version if available, else use raster version for vector outputs
    try:
        vector_version = str(argv[3])    
    except:
        vector_version = raster_version
    
    #####################################
    # Connect to datacube, Dask cluster #
    #####################################
   
    # Connect to datacube    
    dc = datacube.Datacube(app='DEACoastlines_generation')
    
    # Start local dask client
    client = start_local_dask(mem_safety_margin='3gb')
    print(client)    

    ###########################
    # Load supplementary data #
    ###########################

    # Tide points are used to model tides across the extent of the satellite data
    points_gdf = gpd.read_file('input_data/tide_points_coastal.geojson')

    # Albers grid cells used to process the analysis
    studyarea_path = 'input_data/50km_albers_grid_clipped.geojson'
    gridcell_gdf = (gpd.read_file(studyarea_path)
                    .to_crs(epsg=4326)
                    .set_index('id')
                    .loc[[study_area]])

    ################
    # Loading data #
    ################
    
    # Create query
    geopoly = Geometry(gridcell_gdf.iloc[0].geometry, crs=gridcell_gdf.crs)
    query = {'geopolygon': geopoly.buffer(0.05),
             'time': ('1987', '2021'),
             'dask_chunks': {'time': 1, 'x': 3000, 'y': 3000}}

    # Load virtual product    
    ds = load_mndwi(dc, 
                    query, 
                    yaml_path='deacoastlines_virtual_products_v1.0.0.yaml',
                    product_name='ls_nbart_mndwi')
    
    # Temporary workaround map_overlap issues by rechunking
    # if smallest chunk is less than 10
    if ((len(ds.x) % 2000) <= 10) or ((len(ds.y) % 2000) <= 10):
        ds = ds.chunk({'time': 1, 'x': 3200, 'y': 3200})

    # Extract boolean mask
    mask = odc.algo.enum_to_bool(ds.fmask, 
                                 categories=['nodata', 'cloud', 'shadow', 'snow'])

    # Close mask to remove small holes in cloud, open mask to 
    # remove narrow false positive cloud, then dilate
    mask = odc.algo.binary_closing(mask, 2)
    mask_cleaned = odc.algo.mask_cleanup(mask, r=(10, 10))

    # Add new mask as nodata pixels
    ds = odc.algo.erase_bad(ds, mask_cleaned, nodata=np.nan)

    ###################
    # Tidal modelling #
    ###################
    
    # Model tides at point locations
    tidepoints_gdf = model_tides(ds, points_gdf)
    
    # Test if there is data and skip rest of the analysis if not
    if tidepoints_gdf.geometry.unique().shape[0] <= 1:
        sys.exit('Gridcell has 1 or less tidal points; cannot interpolate data')

    # Interpolate tides for each timestep into the spatial extent of the data 
    pool = multiprocessing.Pool(multiprocessing.cpu_count() - 1)
    print(f'Parallelising {multiprocessing.cpu_count() - 1} processes')
    out_list = pool.map(partial(interpolate_tide,
                                tidepoints_gdf=tidepoints_gdf,
                                factor=50), 
                        iterable=[(group.x.values, 
                                   group.y.values, 
                                   group.time.values) 
                                  for (i, group) in ds.groupby('time')])

    # Combine to match the original dataset
    ds['tide_m'] = xr.concat(out_list, dim=ds['time'])    

    # Determine tide cutoff
    tide_cutoff_buff = (
        (ds['tide_m'].max(dim='time') - ds['tide_m'].min(dim='time')) * 0.25)
    tide_cutoff_min = 0.0 - tide_cutoff_buff
    tide_cutoff_max = 0.0 + tide_cutoff_buff
    
    ##############################
    # Generate yearly composites #
    ##############################
    
    # If output folder doesn't exist, create it
    output_dir = f'output_data/{study_area}_{raster_version}'
    os.makedirs(output_dir, exist_ok=True)

    # Iterate through each year and export annual and 3-year gapfill composites
    export_annual_gapfill(ds, 
                          output_dir, 
                          tide_cutoff_min, 
                          tide_cutoff_max)    

    print(f'{(datetime.datetime.now() - start_time).seconds / 60:.1f} minutes')
    
    ##################
    # Run statistics #
    ##################
    
    # Once all rasters have been generated, compute contours and statistics
    os.system(f'python /g/data/r78/DEACoastlines/deacoastlines_statistics.py {study_area} {raster_version} {vector_version}')
    
        
if __name__ == "__main__":
    main()
