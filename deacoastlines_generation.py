#!/usr/bin/env python
# coding: utf-8

import os
import sys
import mock
import otps
import datacube
import datetime
import multiprocessing
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import scipy.interpolate
from affine import Affine
from functools import partial
from shapely.geometry import shape
from datacube.utils.cog import write_cog
from datacube.utils.dask import start_local_dask
from datacube.utils.geometry import GeoBox, Geometry, CRS
from datacube.virtual import catalog_from_file, construct

from collections import Counter
import odc.algo

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


start_time = datetime.datetime.now()



def mostcommon_crs(dc, product, query):
    """
    Takes a given query and returns the most common CRS for observations
    returned for that spatial extent. This can be useful when your study
    area lies on the boundary of two UTM zones, forcing you to decide
    which CRS to use for your `output_crs` in `dc.load`.

    Parameters
    ----------
    dc : datacube Datacube object
        The Datacube to connect to, i.e. `dc = datacube.Datacube()`.
        This allows you to also use development datacubes if required.
    product : str
        A product name to load CRSs from
    query : dict
        A datacube query including x, y and time range to assess for the
        most common CRS

    Returns
    -------
    A EPSG string giving the most common CRS from all datasets returned
    by the query above

    """

    # List of matching products
    matching_datasets = dc.find_datasets(product=product, **query)

    # Extract all CRSs
    crs_list = [str(i.crs) for i in matching_datasets]

    # Identify most common CRS
    crs_counts = Counter(crs_list)
    crs_mostcommon = crs_counts.most_common(1)[0][0]

    # Warn user if multiple CRSs are encountered
    if len(crs_counts.keys()) > 1:

        warnings.warn(f'Multiple UTM zones {list(crs_counts.keys())} '
                      f'were returned for this query. Defaulting to '
                      f'the most common zone: {crs_mostcommon}',
                      UserWarning)

    return crs_mostcommon


def load_ard(dc,
             products=None,
             min_gooddata=0.0,
             fmask_categories=['valid', 'snow', 'water'],
             mask_pixel_quality=True,
             mask_contiguity=False,
             ls7_slc_off=True,
             predicate=None,
             dtype='auto',
             **kwargs):

    """
    Loads and combines Landsat Collection 3 or Sentinel 2 Definitive
    and Near Real Time data for multiple sensors (i.e. ls5t, ls7e and
    ls8c for Landsat; s2a and s2b for Sentinel 2), optionally applies
    pixel quality and contiguity masks, and drops time steps that
    contain greater than a minimum proportion of good quality (e.g. non-
    cloudy or shadowed) pixels.

    The function supports loading the following DEA products:

        ga_ls5t_ard_3
        ga_ls7e_ard_3
        ga_ls8c_ard_3
        s2a_ard_granule
        s2b_ard_granule
        s2a_nrt_granule
        s2b_nrt_granule

    Last modified: June 2020

    Parameters
    ----------
    dc : datacube Datacube object
        The Datacube to connect to, i.e. `dc = datacube.Datacube()`.
        This allows you to also use development datacubes if required.
    products : list
        A list of product names to load data from. Valid options are
        ['ga_ls5t_ard_3', 'ga_ls7e_ard_3', 'ga_ls8c_ard_3'] for Landsat,
        ['s2a_ard_granule', 's2b_ard_granule'] for Sentinel 2 Definitive,
        and ['s2a_nrt_granule', 's2b_nrt_granule'] for Sentinel 2 Near
        Real Time (on the DEA Sandbox only).
    min_gooddata : float, optional
        An optional float giving the minimum percentage of good quality
        pixels required for a satellite observation to be loaded.
        Defaults to 0.0 which will return all observations regardless of
        pixel quality (set to e.g. 0.99 to return only observations with
        more than 99% good quality pixels).
    fmask_categories : list, optional
        An optional list of fmask category names to treat as good
        quality pixels in the above `min_gooddata` calculation, and for
        masking data by pixel quality (if `mask_pixel_quality=True`).
        The default is `['valid', 'snow', 'water']` which will return
        non-cloudy or shadowed land, snow and water pixels. Choose from:
        'nodata', 'valid', 'cloud', 'shadow', 'snow', and 'water'.
    mask_pixel_quality : bool, optional
        An optional boolean indicating whether to mask out poor quality
        pixels using fmask based on the `fmask_categories` provided
        above. The default is True, which will set poor quality pixels
        to NaN if `dtype='auto'` (which will convert the data to
        'float32'), or set poor quality pixels to the data's native
        nodata value if `dtype='native' (which can be useful for
        reducing memory).
    mask_contiguity : str or bool, optional
        An optional string or boolean indicating whether to mask out
        pixels missing data in any band (i.e. "non-contiguous" values).
        This can be important for generating clean composite datasets.
        The default is False, which will ignore non-contiguous values
        completely. If loading NBART data, set the parameter to:
        `mask_contiguity='nbart_contiguity'`. If loading NBAR data,
        specify `mask_contiguity='nbar_contiguity'` instead.
        Non-contiguous pixels will be set to NaN if `dtype='auto'`, or
        set to the data's native nodata value if `dtype='native'`
        (which can be useful for reducing memory).
    dtype : string, optional
        An optional parameter that controls the data type/dtype that
        layers are coerced to after loading. Valid values: 'native',
        'auto', 'float{16|32|64}'. When 'auto' is used, the data will be
        converted to `float32` if masking is used, otherwise data will
        be returned in the native data type of the data. Be aware that
        if data is loaded in its native dtype, nodata and masked
        pixels will be returned with the data's native nodata value
        (typically -999), not NaN.
    ls7_slc_off : bool, optional
        An optional boolean indicating whether to include data from
        after the Landsat 7 SLC failure (i.e. SLC-off). Defaults to
        True, which keeps all Landsat 7 observations > May 31 2003.
    predicate : function, optional
        An optional function that can be passed in to restrict the
        datasets that are loaded by the function. A predicate function
        should take a `datacube.model.Dataset` object as an input (i.e.
        as returned from `dc.find_datasets`), and return a boolean.
        For example, a predicate function could be used to return True
        for only datasets acquired in January:
        `dataset.time.begin.month == 1`
    **kwargs :
        A set of keyword arguments to `dc.load` that define the
        spatiotemporal query and load parameters used to extract data.
        Keyword arguments can either be listed directly in the
        `load_ard` call like any other parameter (e.g.
        `measurements=['nbart_red']`), or by passing in a query kwarg
        dictionary (e.g. `**query`). Keywords can include `measurements`,
        `x`, `y`, `time`, `resolution`, `resampling`, `group_by`, `crs`;
        see the `dc.load` documentation for all possible options:
        https://datacube-core.readthedocs.io/en/latest/dev/api/generate/datacube.Datacube.load.html

    Returns
    -------
    combined_ds : xarray Dataset
        An xarray dataset containing only satellite observations that
        contains greater than `min_gooddata` proportion of good quality
        pixels.

    """
    
    def _dc_query_only(**kw):
        """
        Remove load-only parameters, the rest can be passed to Query

        Returns
        -------
        dict of query parameters
        """

        def _impl(measurements=None,
                  output_crs=None,
                  resolution=None,
                  resampling=None,
                  skip_broken_datasets=None,
                  dask_chunks=None,
                  fuse_func=None,
                  align=None,
                  datasets=None,
                  progress_cbk=None,
                  group_by=None,
                  **query):
            return query

        return _impl(**kw)


    def _common_bands(dc, products):
        """
        Takes a list of products and returns a list of measurements/bands
        that are present in all products

        Returns
        -------
        List of band names
        """
        common = None
        bands = None

        for p in products:
            p = dc.index.products.get_by_name(p)
            if common is None:
                common = set(p.measurements)
                bands = list(p.measurements)
            else:
                common = common.intersection(set(p.measurements))
        return [band for band in bands if band in common]


    def _dc_query_only(**kw):
        """
        Remove load-only parameters, the rest can be passed to Query

        Returns
        -------
        dict of query parameters
        """

        def _impl(measurements=None,
                  output_crs=None,
                  resolution=None,
                  resampling=None,
                  skip_broken_datasets=None,
                  dask_chunks=None,
                  fuse_func=None,
                  align=None,
                  datasets=None,
                  progress_cbk=None,
                  group_by=None,
                  **query):
            return query

        return _impl(**kw)


    def _common_bands(dc, products):
        """
        Takes a list of products and returns a list of measurements/bands
        that are present in all products

        Returns
        -------
        List of band names
        """
        common = None
        bands = None

        for p in products:
            p = dc.index.products.get_by_name(p)
            if common is None:
                common = set(p.measurements)
                bands = list(p.measurements)
            else:
                common = common.intersection(set(p.measurements))
        return [band for band in bands if band in common]

    #########
    # Setup #
    #########

    # Use 'nbart_contiguity' by default if mask_contiguity is true
    if mask_contiguity is True:
        mask_contiguity = 'nbart_contiguity'

    # We deal with `dask_chunks` separately
    dask_chunks = kwargs.pop('dask_chunks', None)
    requested_measurements = kwargs.pop('measurements', None)

    # Warn user if they combine lazy load with min_gooddata
    if (min_gooddata > 0.0) and dask_chunks is not None:
        warnings.warn("Setting 'min_gooddata' percentage to > 0.0 "
                      "will cause dask arrays to compute when "
                      "loading pixel-quality data to calculate "
                      "'good pixel' percentage. This can "
                      "slow the return of your dataset.")

    # Verify that products were provided, and determine if Sentinel-2
    # or Landsat data is being loaded
    if not products:
        raise ValueError("Please provide a list of product names "
                         "to load data from. Valid options are: \n"
                         "['ga_ls5t_ard_3', 'ga_ls7e_ard_3', 'ga_ls8c_ard_3'] "
                         "for Landsat, ['s2a_ard_granule', "
                         "'s2b_ard_granule'] \nfor Sentinel 2 Definitive, or "
                         "['s2a_nrt_granule', 's2b_nrt_granule'] for "
                         "Sentinel 2 Near Real Time")
    elif all(['ls' in product for product in products]):
        product_type = 'ls'
    elif all(['s2' in product for product in products]):
        product_type = 's2'

    fmask_band = 'fmask'
    measurements = (requested_measurements.copy() if
                    requested_measurements else None)

    if measurements is None:

        # Deal with "load all" case: pick a set of bands common across
        # all products
        measurements = _common_bands(dc, products)

        # If no `measurements` are specified, Landsat ancillary bands are
        # loaded with a 'oa_' prefix, but Sentinel-2 bands are not. As a
        # work-around, we need to rename the default contiguity and fmask
        # bands if loading Landsat data without specifying `measurements`
        if product_type == 'ls':
            mask_contiguity = (f'oa_{mask_contiguity}' if
                               mask_contiguity else False)
            fmask_band = f'oa_{fmask_band}'

    # If `measurements` are specified but do not include fmask or
    # contiguity variables, add these to `measurements`
    if fmask_band not in measurements:
        measurements.append(fmask_band)
    if mask_contiguity and mask_contiguity not in measurements:
        measurements.append(mask_contiguity)

    # Get list of data and mask bands so that we can later exclude
    # mask bands from being masked themselves
    data_bands = [band for band in measurements if
                  band not in (fmask_band, mask_contiguity)]
    mask_bands = [band for band in measurements if band not in data_bands]

    #################
    # Find datasets #
    #################

    # Pull out query params only to pass to dc.find_datasets
    query = _dc_query_only(**kwargs)

    # Extract datasets for each product using subset of dcload_kwargs
    dataset_list = []

    # Get list of datasets for each product
    print('Finding datasets')
    for product in products:

        # Obtain list of datasets for product
        print(f'    {product} (ignoring SLC-off observations)'
              if not ls7_slc_off and product == 'ga_ls7e_ard_3'
              else f'    {product}')
        datasets = dc.find_datasets(product=product, **query)

        # Remove Landsat 7 SLC-off observations if ls7_slc_off=False
        if not ls7_slc_off and product == 'ga_ls7e_ard_3':
            datasets = [i for i in datasets if
                        normalise_dt(i.time.begin) <
                        datetime.datetime(2003, 5, 31)]

        # Add any returned datasets to list
        dataset_list.extend(datasets)

    # Raise exception if no datasets are returned
    if len(dataset_list) == 0:
        raise ValueError("No data available for query: ensure that "
                         "the products specified have data for the "
                         "time and location requested")

    # If predicate is specified, use this function to filter the list
    # of datasets prior to load
    if predicate:
        print(f'Filtering datasets using predicate function')
        dataset_list = [ds for ds in dataset_list if predicate(ds)]

    # Raise exception if filtering removes all datasets
    if len(dataset_list) == 0:
        raise ValueError("No data available after filtering with "
                         "predicate function")

    #############
    # Load data #
    #############

    # Note we always load using dask here so that we can lazy load data
    # before filtering by good data
    ds = dc.load(datasets=dataset_list,
                 measurements=measurements,
                 dask_chunks={} if dask_chunks is None else dask_chunks,
                 **kwargs)

    ####################
    # Filter good data #
    ####################

    # Calculate pixel quality mask
    pq_mask = odc.algo.fmask_to_bool(ds[fmask_band],
                                     categories=fmask_categories)

    # The good data percentage calculation has to load in all `fmask`
    # data, which can be slow. If the user has chosen no filtering
    # by using the default `min_gooddata = 0`, we can skip this step
    # completely to save processing time
    if min_gooddata > 0.0:

        # Compute good data for each observation as % of total pixels
        print('Counting good quality pixels for each time step')
        data_perc = (pq_mask.sum(axis=[1, 2], dtype='int32') /
                     (pq_mask.shape[1] * pq_mask.shape[2]))
        keep = data_perc >= min_gooddata

        # Filter by `min_gooddata` to drop low quality observations
        total_obs = len(ds.time)
        ds = ds.sel(time=keep)
        pq_mask = pq_mask.sel(time=keep)

        print(f'Filtering to {len(ds.time)} out of {total_obs} '
              f'time steps with at least {min_gooddata:.1%} '
              f'good quality pixels')

    ###############
    # Apply masks #
    ###############

    # Create an overall mask to hold both pixel quality and contiguity
    mask = None

    # Add pixel quality mask to overall mask
    if mask_pixel_quality:
        print('Applying pixel quality/cloud mask')
        mask = pq_mask

    # Add contiguity mask to overall mask
    if mask_contiguity:
        print('Applying contiguity mask')
        cont_mask = ds[mask_contiguity] == 1

        # If mask already has data if mask_pixel_quality == True,
        # multiply with cont_mask to perform a logical 'or' operation
        # (keeping only pixels good in both)
        mask = cont_mask if mask is None else mask * cont_mask

    # Split into data/masks bands, as conversion to float and masking
    # should only be applied to data bands
    ds_data = ds[data_bands]
    ds_masks = ds[mask_bands]

    # Mask data if either of the above masks were generated
    if mask is not None:
        ds_data = odc.algo.keep_good_only(ds_data, where=mask)

    # Automatically set dtype to either native or float32 depending
    # on whether masking was requested
    if dtype == 'auto':
        dtype = 'native' if mask is None else 'float32'

    # Set nodata values using odc.algo tools to reduce peak memory
    # use when converting data dtype
    if dtype != 'native':
        ds_data = odc.algo.to_float(ds_data, dtype=dtype)

    # Put data and mask bands back together
    attrs = ds.attrs
    ds = xr.merge([ds_data, ds_masks])
    ds.attrs.update(attrs)

    ###############
    # Return data #
    ###############

    # Drop bands not originally requested by user
    if requested_measurements:
        ds = ds[requested_measurements]

    # If user supplied dask_chunks, return data as a dask array without
    # actually loading it in
    if dask_chunks is not None:
        print(f'Returning {len(ds.time)} time steps as a dask array')
        return ds
    else:
        print(f'Loading {len(ds.time)} time steps')
        return ds.compute()
    
    
def calculate_indices(ds,
                      index=None,
                      collection=None,
                      custom_varname=None,
                      normalise=True,
                      drop=False,
                      deep_copy=True):
    """
    Takes an xarray dataset containing spectral bands, calculates one of
    a set of remote sensing indices, and adds the resulting array as a 
    new variable in the original dataset.  
    
    Last modified: September 2019
    
    Parameters
    ----------  
    ds : xarray Dataset
        A two-dimensional or multi-dimensional array with containing the 
        spectral bands required to calculate the index. These bands are 
        used as inputs to calculate the selected water index.
    index : str or list of strs
        A string giving the name of the index to calculate or a list of 
        strings giving the names of the indices to calculate:
        'AWEI_ns (Automated Water Extraction Index,
                  no shadows, Feyisa 2014)
        'AWEI_sh' (Automated Water Extraction Index,
                   shadows, Feyisa 2014)
        'BAEI' (Built-Up Area Extraction Index, Bouzekri et al. 2015) 
        'BAI' (Burn Area Index, Martin 1998)
        'BSI' (Bare Soil Index, Rikimaru et al. 2002)
        'BUI' (Built-Up Index, He et al. 2010)
        'CMR' (Clay Minerals Ratio, Drury 1987)
        'EVI' (Enhanced Vegetation Index, Huete 2002)
        'FMR' (Ferrous Minerals Ratio, Segal 1982)
        'IOR' (Iron Oxide Ratio, Segal 1982)  
        'LAI' (Leaf Area Index, Boegh 2002)
        'MNDWI' (Modified Normalised Difference Water Index, Xu 1996) 
        'MSAVI' (Modified Soil Adjusted Vegetation Index, 
                 Qi et al. 1994)              
        'NBI' (New Built-Up Index, Jieli et al. 2010)
        'NBR' (Normalised Burn Ratio, Lopez Garcia 1991)
        'NDBI' (Normalised Difference Built-Up Index, Zha 2003)
        'NDCI' (Normalised Difference Chlorophyll Index, 
                Mishra & Mishra, 2012)
        'NDMI' (Normalised Difference Moisture Index, Gao 1996)        
        'NDSI' (Normalised Difference Snow Index, Hall 1995)
        'NDVI' (Normalised Difference Vegetation Index, Rouse 1973)
        'NDWI' (Normalised Difference Water Index, McFeeters 1996)
        'SAVI' (Soil Adjusted Vegetation Index, Huete 1988)
        'TCB' (Tasseled Cap Brightness, Crist 1985)
        'TCG' (Tasseled Cap Greeness, Crist 1985)
        'TCW' (Tasseled Cap Wetness, Crist 1985)
        'WI' (Water Index, Fisher 2016) 
    collection : str
        An string that tells the function what data collection is 
        being used to calculate the index. This is necessary because 
        different collections use different names for bands covering 
        a similar spectra. Valid options are 'ga_ls_2' (for GA 
        Landsat Collection 2), 'ga_ls_3' (for GA Landsat Collection 3) 
        and 'ga_s2_1' (for GA Sentinel 2 Collection 1).
    custom_varname : str, optional
        By default, the original dataset will be returned with 
        a new index variable named after `index` (e.g. 'NDVI'). To 
        specify a custom name instead, you can supply e.g. 
        `custom_varname='custom_name'`. Defaults to None, which uses
        `index` to name the variable. 
    normalise : bool, optional
        Some coefficient-based indices (e.g. 'WI', 'BAEI', 'AWEI_ns', 
        'AWEI_sh', 'TCW', 'TCG', 'TCB', 'EVI', 'LAI', 'SAVI', 'MSAVI') 
        produce different results if surface reflectance values are not 
        scaled between 0.0 and 1.0 prior to calculating the index. 
        Setting `normalise=True` first scales values to a 0.0-1.0 range
        by dividing by 10000.0. Defaults to True.  
    drop : bool, optional
        Provides the option to drop the original input data, thus saving 
        space. if drop = True, returns only the index and its values.
    deep_copy: bool, optional
        If deep_copy=False, calculate_indices will modify the original
        array, adding bands to the input dataset and not removing them.
        If the calculate_indices function is run more than once, variables
        may be dropped incorrectly producing unexpected behaviour. This is
        a bug and may be fixed in future releases. This is only a problem 
        when drop=True.
    
        
    Returns
    -------
    ds : xarray Dataset
        The original xarray Dataset inputted into the function, with a 
        new varible containing the remote sensing index as a DataArray.
        If drop = True, the new variable/s as DataArrays in the 
        original Dataset. 
    """
    
    # Set ds equal to a copy of itself in order to prevent the function 
    # from editing the input dataset. This is to prevent unexpected 
    # behaviour though it uses twice as much memory.    
    if deep_copy:
        ds = ds.copy(deep=True)
    
    # Capture input band names in order to drop these if drop=True
    if drop:
        bands_to_drop=list(ds.data_vars)
        print(f'Dropping bands {bands_to_drop}')

    # Dictionary containing remote sensing index band recipes
    index_dict = {
                  # Normalised Difference Vegation Index, Rouse 1973
                  'NDVI': lambda ds: (ds.nir - ds.red) /
                                     (ds.nir + ds.red),

                  # Enhanced Vegetation Index, Huete 2002
                  'EVI': lambda ds: ((2.5 * (ds.nir - ds.red)) /
                                     (ds.nir + 6 * ds.red -
                                      7.5 * ds.blue + 1)),

                  # Leaf Area Index, Boegh 2002
                  'LAI': lambda ds: (3.618 * ((2.5 * (ds.nir - ds.red)) /
                                     (ds.nir + 6 * ds.red -
                                      7.5 * ds.blue + 1)) - 0.118),

                  # Soil Adjusted Vegetation Index, Huete 1988
                  'SAVI': lambda ds: ((1.5 * (ds.nir - ds.red)) /
                                      (ds.nir + ds.red + 0.5)),
      
                  # Mod. Soil Adjusted Vegetation Index, Qi et al. 1994
                  'MSAVI': lambda ds: ((2 * ds.nir + 1 - 
                                      ((2 * ds.nir + 1)**2 - 
                                       8 * (ds.nir - ds.red))**0.5) / 2),    

                  # Normalised Difference Moisture Index, Gao 1996
                  'NDMI': lambda ds: (ds.nir - ds.swir1) /
                                     (ds.nir + ds.swir1),

                  # Normalised Burn Ratio, Lopez Garcia 1991
                  'NBR': lambda ds: (ds.nir - ds.swir2) /
                                    (ds.nir + ds.swir2),

                  # Burn Area Index, Martin 1998
                  'BAI': lambda ds: (1.0 / ((0.10 - ds.red) ** 2 +
                                            (0.06 - ds.nir) ** 2)),
        
                 # Normalised Difference Chlorophyll Index, 
                 # (Mishra & Mishra, 2012)
                  'NDCI': lambda ds: (ds.red_edge_1 - ds.red) /
                                     (ds.red_edge_1 + ds.red),

                  # Normalised Difference Snow Index, Hall 1995
                  'NDSI': lambda ds: (ds.green - ds.swir1) /
                                     (ds.green + ds.swir1),

                  # Normalised Difference Water Index, McFeeters 1996
                  'NDWI': lambda ds: (ds.green - ds.nir) /
                                     (ds.green + ds.nir),

                  # Modified Normalised Difference Water Index, Xu 2006
                  'MNDWI': lambda ds: (ds.green - ds.swir1) /
                                      (ds.green + ds.swir1),
      
                  # Normalised Difference Built-Up Index, Zha 2003
                  'NDBI': lambda ds: (ds.swir1 - ds.nir) /
                                     (ds.swir1 + ds.nir),
      
                  # Built-Up Index, He et al. 2010
                  'BUI': lambda ds:  ((ds.swir1 - ds.nir) /
                                      (ds.swir1 + ds.nir)) -
                                     ((ds.nir - ds.red) /
                                      (ds.nir + ds.red)),
      
                  # Built-up Area Extraction Index, Bouzekri et al. 2015
                  'BAEI': lambda ds: (ds.red + 0.3) /
                                     (ds.green + ds.swir1),
      
                  # New Built-up Index, Jieli et al. 2010
                  'NBI': lambda ds: (ds.swir1 + ds.red) / ds.nir,
      
                  # Bare Soil Index, Rikimaru et al. 2002
                  'BSI': lambda ds: ((ds.swir1 + ds.red) - 
                                     (ds.nir + ds.blue)) / 
                                    ((ds.swir1 + ds.red) + 
                                     (ds.nir + ds.blue)),

                  # Automated Water Extraction Index (no shadows), Feyisa 2014
                  'AWEI_ns': lambda ds: (4 * (ds.green - ds.swir1) -
                                        (0.25 * ds.nir * + 2.75 * ds.swir2)),

                  # Automated Water Extraction Index (shadows), Feyisa 2014
                  'AWEI_sh': lambda ds: (ds.blue + 2.5 * ds.green -
                                         1.5 * (ds.nir + ds.swir1) -
                                         0.25 * ds.swir2),

                  # Water Index, Fisher 2016
                  'WI': lambda ds: (1.7204 + 171 * ds.green + 3 * ds.red -
                                    70 * ds.nir - 45 * ds.swir1 -
                                    71 * ds.swir2),

                  # Tasseled Cap Wetness, Crist 1985
                  'TCW': lambda ds: (0.0315 * ds.blue + 0.2021 * ds.green +
                                     0.3102 * ds.red + 0.1594 * ds.nir +
                                    -0.6806 * ds.swir1 + -0.6109 * ds.swir2),

                  # Tasseled Cap Greeness, Crist 1985
                  'TCG': lambda ds: (-0.1603 * ds.blue + -0.2819 * ds.green +
                                     -0.4934 * ds.red + 0.7940 * ds.nir +
                                     -0.0002 * ds.swir1 + -0.1446 * ds.swir2),

                  # Tasseled Cap Brightness, Crist 1985
                  'TCB': lambda ds: (0.2043 * ds.blue + 0.4158 * ds.green +
                                     0.5524 * ds.red + 0.5741 * ds.nir +
                                     0.3124 * ds.swir1 + -0.2303 * ds.swir2),

                  # Clay Minerals Ratio, Drury 1987
                  'CMR': lambda ds: (ds.swir1 / ds.swir2),

                  # Ferrous Minerals Ratio, Segal 1982
                  'FMR': lambda ds: (ds.swir1 / ds.nir),

                  # Iron Oxide Ratio, Segal 1982
                  'IOR': lambda ds: (ds.red / ds.blue)
    }
    
    # If index supplied is not a list, convert to list. This allows us to
    # iterate through either multiple or single indices in the loop below
    indices = index if isinstance(index, list) else [index]
    
    #calculate for each index in the list of indices supplied (indexes)
    for index in indices:

        # Select an index function from the dictionary
        index_func = index_dict.get(str(index))

        # If no index is provided or if no function is returned due to an 
        # invalid option being provided, raise an exception informing user to 
        # choose from the list of valid options
        if index is None:

            raise ValueError(f"No remote sensing `index` was provided. Please "
                              "refer to the function \ndocumentation for a full "
                              "list of valid options for `index` (e.g. 'NDVI')")

        elif (index in ['WI', 'BAEI', 'AWEI_ns', 'AWEI_sh', 'TCW', 
                        'TCG', 'TCB', 'EVI', 'LAI', 'SAVI', 'MSAVI'] 
              and not normalise):

            warnings.warn(f"\nA coefficient-based index ('{index}') normally "
                           "applied to surface reflectance values in the \n"
                           "0.0-1.0 range was applied to values in the 0-10000 "
                           "range. This can produce unexpected results; \nif "
                           "required, resolve this by setting `normalise=True`")

        elif index_func is None:

            raise ValueError(f"The selected index '{index}' is not one of the "
                              "valid remote sensing index options. \nPlease "
                              "refer to the function documentation for a full "
                              "list of valid options for `index`")

        # Rename bands to a consistent format if depending on what collection
        # is specified in `collection`. This allows the same index calculations
        # to be applied to all collections. If no collection was provided, 
        # raise an exception.
        if collection is None:

            raise ValueError("'No `collection` was provided. Please specify "
                             "either 'ga_ls_2', 'ga_ls_3' or 'ga_s2_1' \nto "
                             "ensure the function calculates indices using the "
                             "correct spectral bands")

        elif collection == 'ga_ls_3':

            # Dictionary mapping full data names to simpler 'red' alias names
            bandnames_dict = {
                'nbart_nir': 'nir',
                'nbart_red': 'red',
                'nbart_green': 'green',
                'nbart_blue': 'blue',
                'nbart_swir_1': 'swir1',
                'nbart_swir_2': 'swir2',
                'nbar_red': 'red',
                'nbar_green': 'green',
                'nbar_blue': 'blue',
                'nbar_nir': 'nir',
                'nbar_swir_1': 'swir1',
                'nbar_swir_2': 'swir2'
            }

            # Rename bands in dataset to use simple names (e.g. 'red')
            bands_to_rename = {
                a: b for a, b in bandnames_dict.items() if a in ds.variables
            }

        elif collection == 'ga_s2_1':

            # Dictionary mapping full data names to simpler 'red' alias names
            bandnames_dict = {
                'nbart_red': 'red',
                'nbart_green': 'green',
                'nbart_blue': 'blue',
                'nbart_nir_1': 'nir',
                'nbart_red_edge_1': 'red_edge_1', 
                'nbart_red_edge_2': 'red_edge_2',    
                'nbart_swir_2': 'swir1',
                'nbart_swir_3': 'swir2',
                'nbar_red': 'red',
                'nbar_green': 'green',
                'nbar_blue': 'blue',
                'nbar_nir_1': 'nir',
                'nbar_red_edge_1': 'red_edge_1',   
                'nbar_red_edge_2': 'red_edge_2',   
                'nbar_swir_2': 'swir1',
                'nbar_swir_3': 'swir2'
            }

            # Rename bands in dataset to use simple names (e.g. 'red')
            bands_to_rename = {
                a: b for a, b in bandnames_dict.items() if a in ds.variables
            }

        elif collection == 'ga_ls_2':

            # Pass an empty dict as no bands need renaming
            bands_to_rename = {}

        # Raise error if no valid collection name is provided:
        else:
            raise ValueError(f"'{collection}' is not a valid option for "
                              "`collection`. Please specify either \n"
                              "'ga_ls_2', 'ga_ls_3' or 'ga_s2_1'")

        # Apply index function 
        try:
            # If normalised=True, divide data by 10,000 before applying func
            mult = 10000.0 if normalise else 1.0
            index_array = index_func(ds.rename(bands_to_rename) / mult)
        except AttributeError:
            raise ValueError(f'Please verify that all bands required to '
                             f'compute {index} are present in `ds`. \n'
                             f'These bands may vary depending on the `collection` '
                             f'(e.g. the Landsat `nbart_nir` band \n'
                             f'is equivelent to `nbart_nir_1` for Sentinel 2)')

        # Add as a new variable in dataset
        output_band_name = custom_varname if custom_varname else index
        ds[output_band_name] = index_array
    
    # Once all indexes are calculated, drop input bands if drop=True
    if drop: 
        ds = ds.drop(bands_to_drop)

    # Return input dataset with added water index variable
    return ds


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
    return GeoBox(width=shape[1], height=shape[0], affine=affine, crs=crs)


def load_mndwi(dc, 
               query,
               yaml_path, 
               product_name='ls_nbart_mndwi',
               virtual_products=True):
    """
    This function uses virtual products to load data from GA Collection 
    3 Landsat 5, 7 and 8, calculate custom remote sensing indices, and 
    return the data as a single xarray.Dataset.
    
    To minimise resampling effects and maintain the highest data 
    fidelity required for subpixel coastline extraction, this workflow 
    applies masking and index calculation at native resolution, and 
    only re-projects to the most common CRS for the query using average 
    resampling in the final step.
    """

    # Identify the most common CRS in the region, so data can be loaded with 
    # minimal distortion. The dictionary comprehension is required as 
    # dc.find_datasets does not work in combination with dask_chnks
    crs = mostcommon_crs(dc=dc, product='ga_ls5t_ard_3', 
                         query={k: v for k, v in query.items() if 
                                k not in ['dask_chunks']})
    
    if virtual_products:
    
        # Load in virtual product catalogue and select MNDWI product
        catalog = catalog_from_file(yaml_path)
        product = catalog[product_name]

        # Construct a new version of the product using most common CRS
        product_reproject = construct(input=product,
                                      reproject={'output_crs': str(crs), 
                                                 'resolution': (-30, 30),
                                                 'align': (15, 15)},          
                                      resampling='average')

        # Determine geobox with custom function to increase lazy loading 
        # speed (will eventually be done automatically within virtual 
        # products)
        with mock.patch('datacube.virtual.impl.native_geobox', 
                        side_effect=custom_native_geobox):
            ds = product_reproject.load(dc, **query)
    
    else:
               
        ds = load_ard(dc=dc, 
              measurements=['nbart_blue', 'nbart_green', 'nbart_red', 
                            'nbart_nir', 'nbart_swir_1', 'nbart_swir_2'], 
              min_gooddata=0.0,
              products=['ga_ls5t_ard_3', 'ga_ls7e_ard_3', 'ga_ls8c_ard_3'], 
              output_crs=crs,
              resampling={'fmask': 'nearest', 
                          'oa_fmask': 'nearest', 
                          'nbart_contiguity': 'nearest',
                          'oa_nbart_contiguity': 'nearest',
                          '*': 'cubic'},
              resolution=(-30, 30),  
              gqa_iterative_mean_xy=[0, 1],
              align=(15, 15),
              group_by='solar_day',
              mask_contiguity=False,
              **query)

        ds = (calculate_indices(ds, index=['MNDWI'], 
                                collection='ga_ls_3', 
                                drop=True)
              .rename({'MNDWI': 'mndwi'}))
        
        
    return ds


def model_tides(ds, points_gdf, extent_buffer=0.05):
    """
    Takes an xarray.Dataset (`ds`), extracts a subset of tide modelling 
    points from a geopandas.GeoDataFrame based on`ds`'s extent, then 
    uses the OTPS tidal model to model tide heights for every point
    at every time step in `ds`.
    
    The output is a geopandas.GeoDataFrame with a "time" index 
    (matching the time steps in `ds`), and a "tide_m" column giving the 
    tide heights at each point location.
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


def interpolate_2d(x_coords, 
                   y_coords, 
                   z_coords, 
                   grid_x_ds,
                   grid_y_ds,
                   method='linear',
                   factor=1,
                   **kwargs):
    
    """
    This function takes points with X, Y and Z coordinates, and 
    interpolates Z-values across the extent of an existing xarray 
    dataset. This can be useful for producing smooth surfaces from point
    data that can be compared directly against satellite data derived 
    from an OpenDataCube query.
    
    Supported interpolation methods include 'linear', 'nearest' and
    'cubic (using `scipy.interpolate.griddata`), and 'rbf' (using 
    `scipy.interpolate.Rbf`).
    
    Last modified: March 2019
    
    Parameters
    ----------  
    x_coords, y_coords : numpy array
        Arrays containing X and Y coordinates for all points (e.g. 
        longitudes and latitudes).
    z_coords : numpy array
        An array containing Z coordinates for all points (e.g. 
        elevations). These are the values you wish to interpolate 
        between.
    method : string, optional
        The method used to interpolate between point values. This string
        is either passed to `scipy.interpolate.griddata` (for 'linear', 
        'nearest' and 'cubic' methods), or used to specify Radial Basis 
        Function interpolation using `scipy.interpolate.Rbf` ('rbf').
        Defaults to 'linear'.
    factor : int, optional
        An optional integer that can be used to subsample the spatial 
        interpolation extent to obtain faster interpolation times, then
        up-sample this array back to the original dimensions of the 
        data as a final step. For example, setting `factor=10` will 
        interpolate data into a grid that has one tenth of the 
        resolution of `ds`. This approach will be significantly faster 
        than interpolating at full resolution, but will potentially 
        produce less accurate or reliable results.
    **kwargs : 
        Optional keyword arguments to pass to either 
        `scipy.interpolate.griddata` (if `method` is 'linear', 'nearest' 
        or 'cubic'), or `scipy.interpolate.Rbf` (is `method` is 'rbf').
      
    Returns
    -------
    interp_2d_array : xarray DataArray
        An xarray DataArray containing with x and y coordinates copied 
        from `ds_array`, and Z-values interpolated from the points data. 
    """    
  
    # Extract xy and elev points
    points_xy = np.vstack([x_coords, y_coords]).T
    
    # Extract x and y coordinates to interpolate into. 
    # If `factor` is greater than 1, the coordinates will be subsampled 
    # for faster run-times. If the last x or y value in the subsampled 
    # grid aren't the same as the last x or y values in the original 
    # full resolution grid, add the final full resolution grid value to 
    # ensure data is interpolated up to the very edge of the array
    if grid_x_ds[::factor][-1] == grid_x_ds[-1]:
        x_grid_coords = grid_x_ds[::factor]
    else:
        x_grid_coords = grid_x_ds[::factor].tolist() + [grid_x_ds[-1]]
        
    if grid_y_ds[::factor][-1] == grid_y_ds[-1]:
        y_grid_coords = grid_y_ds[::factor]
    else:
        y_grid_coords = grid_y_ds[::factor].tolist() + [grid_y_ds[-1]]

    # Create grid to interpolate into
    grid_y, grid_x = np.meshgrid(x_grid_coords, y_grid_coords)
        
    # Apply scipy.interpolate.griddata interpolation methods
    if method in ('linear', 'nearest', 'cubic'):       

        # Interpolate x, y and z values 
        interp_2d = scipy.interpolate.griddata(points=points_xy, 
                                                values=z_coords, 
                                                xi=(grid_y, grid_x), 
                                                method=method,
                                                **kwargs)
        
    # Apply Radial Basis Function interpolation
    elif method == 'rbf':
        
        # Interpolate x, y and z values 
        rbf = scipy.interpolate.Rbf(x_coords, y_coords, z_coords, **kwargs)  
        interp_2d = rbf(grid_y, grid_x)

    # Create xarray dataarray from the data and resample to ds coords
    interp_2d_da = xr.DataArray(interp_2d,
                                coords=[y_grid_coords, x_grid_coords], 
                                dims=['y', 'x'])
    
    # If factor is greater than 1, resample the interpolated array to
    # match the input `ds` array
    ds_to_interp = xr.DataArray(np.ones(shape=(len(grid_y_ds), len(grid_x_ds))),
             coords=[grid_y_ds, grid_x_ds], 
             dims=['y', 'x'])
    
    if factor > 1: 
        interp_2d_da = interp_2d_da.interp_like(ds_to_interp)

    return interp_2d_da


def interpolate_tide(timestep_tuple, 
                     tidepoints_gdf, 
                     method='rbf', 
                     factor=20):    
    """
    Extract a subset of tide modelling point data for a given time-step,
    then interpolate these tides into the extent of the xarray dataset.
    """  
  
    # Extract subset of observations based on timestamp of imagery
    time_string = str(timestep_tuple[2])[0:19].replace('T', ' ')
    tidepoints_subset = tidepoints_gdf.loc[time_string]
    print(f'{time_string:<80}', end='\r')
    
    # Get lists of x, y and z (tide height) data to interpolate
    x_coords = tidepoints_subset.geometry.x.values.astype('float32')
    y_coords = tidepoints_subset.geometry.y.values.astype('float32')
    z_coords = tidepoints_subset.tide_m.values.astype('float32')    
      
    # Interpolate tides into the extent of the satellite timestep
    out_tide = interpolate_2d(x_coords=x_coords,
                              y_coords=y_coords,
                              z_coords=z_coords,
                              grid_x_ds=timestep_tuple[0],
                              grid_y_ds=timestep_tuple[1],
                              method=method,
                              factor=factor)
    
    # Return data as a Float32 to conserve memory
    return out_tide.astype('float32')


def multiprocess_apply(ds, dim, func):
    """
    Applies a custom function along the dimension of an xarray.Dataset,
    then combines the output to match the original dataset.
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
    """
    
    # Print status
    year = year_ds.time[0].dt.year.item()
    print(f'Processing {year}')
    
    # Determine what pixels were acquired in selected tide range, and 
    # drop time-steps without any relevant pixels to reduce data to load
    tide_bool = ((year_ds.tide_m >= tide_cutoff_min) & 
                 (year_ds.tide_m <= tide_cutoff_max))
    year_ds = year_ds.sel(time=tide_bool.sum(dim=['x', 'y']) > 0)
    
    # Apply mask, and load in corresponding high tide data
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

        str_usage = "You must specify a study area ID and name"
        print(str_usage)
        sys.exit()
        
    # Set study area and name for analysis
    study_area = int(argv[1])
    output_name = str(argv[2])    
   
    # Connect to datacube    
    dc = datacube.Datacube(app='DEACoastLines_generation')
    
    # Start local dask client
    client = start_local_dask(mem_safety_margin='3gb')
    print(client)    

    ###########################
    # Load supplementary data #
    ###########################

    # Tide points are used to model tides across the extent of the satellite data
    points_gdf = gpd.read_file('input_data/tide_points_coastal.geojson')

    # Albers grid cells used to process the analysis
    gridcell_gdf = (gpd.read_file('input_data/50km_albers_grid_clipped.geojson')
                    .to_crs(epsg=4326)
                    .set_index('id')
                    .loc[[study_area]])

    ################
    # Loading data #
    ################
    
    # Create query
    geopoly = Geometry(gridcell_gdf.iloc[0].geometry, crs=gridcell_gdf.crs)
    query = {'geopolygon': geopoly.buffer(0.05),
             'time': ('1987', '2020'),
             'cloud_cover': [0, 90],
             'dask_chunks': {'time': 1, 'x': 2000, 'y': 2000}}

    # Load virtual product    
    ds = load_mndwi(dc, 
                    query, 
                    yaml_path='deacoastlines_virtual_products.yaml',
                    virtual_products=False)
    
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
    output_dir = f'output_data/{study_area}_{output_name}'
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
    os.system(f'python /g/data/r78/DEACoastLines/deacoastlines_statistics.py {study_area} {output_name}')
    
        
if __name__ == "__main__":
    main()