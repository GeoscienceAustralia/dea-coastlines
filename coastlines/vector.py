#!/usr/bin/env python
# coding: utf-8

# This code conducts vector subpixel shoreline extraction for DEA
# Coastlines:
#
#     * Apply morphological extraction algorithms to mask annual median
#       composite rasters to a valid coastal region
#     * Extract waterline vectors using subpixel waterline extraction
#       (Bishop-Taylor et al. 2019b; https://doi.org/10.3390/rs11242984)
#     * Compute rates of coastal change at every 30 m of coastline
#       using linear regression

import glob
import os
import sys
import warnings

import click
import pyproj
import odc.algo
import odc.geo.xr
import numpy as np
import pandas as pd
import xarray as xr
import geohash as gh
import geopandas as gpd
from affine import Affine
from rasterio.features import sieve
from scipy.stats import circstd, circmean, linregress
from shapely.geometry import box
from shapely.ops import nearest_points
from skimage.measure import label, regionprops
from skimage.morphology import (
    black_tophat,
    binary_closing,
    binary_dilation,
    binary_erosion,
    dilation,
    disk,
)

import datacube
from datacube.utils.aws import configure_s3_access

from coastlines.utils import configure_logging, load_config
from dea_tools.spatial import subpixel_contours, xr_vectorize, xr_rasterize

# Hide specific warnings
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=DeprecationWarning)
warnings.simplefilter(action="ignore", category=RuntimeWarning)
pd.options.mode.chained_assignment = None


def load_rasters(
    path,
    raster_version,
    study_area,
    water_index="mndwi",
    start_year=1988,
    end_year=2021,
):
    """
    Loads DEA Coastlines water index (e.g. 'MNDWI'), 'count',
    and 'stdev' rasters for both annual and three-year gapfill data
    into a consistent `xarray.Dataset` format for further analysis.

    Parameters:
    -----------
    path : string
        A string giving the directory containing raster outputs.
    raster_version : string
        A string giving the unique analysis name (e.g. 'v0.3.0') used
        to load raster files.
    study_area : string or int
        A string giving the study area grid cell used to name raster files
        (e.g. tile `6931`).
    water_index : string, optional
        A string giving the name of the water index to load. Defaults
        to 'mndwi', which will load raster files produced using the
        Modified Normalised Difference Water Index.
    start_year : integer, optional
        The first annual layer to include in the analysis. Defaults to
        1988.
    end_year : integer, optional
        The final annual layer to include in the analysis. Defaults to
        2021.

    Returns:
    --------
    yearly_ds : xarray.Dataset
        An `xarray.Dataset` containing annual input rasters.
        The dataset contains water index (e.g. 'MNDWI'),
        'count', and 'stdev' arrays for each year from 1988 onward.
    gapfill_ds : xarray.Dataset
        An `xarray.Dataset` containing three-year gapfill rasters.
        The dataset contains water index (e.g. 'MNDWI'),
        'count', and 'stdev' arrays for each year from
        `start_year` to `end_year`.

    """

    # List to hold output Datasets
    ds_list = []

    for layer_type in [".tif", "_gapfill.tif"]:

        # List to hold output DataArrays
        da_list = []

        for layer_name in [f"{water_index}", "count", "stdev"]:

            # Get paths of files that match pattern
            paths = glob.glob(
                f"{path}/{raster_version}/"
                f"{study_area}_{raster_version}/"
                f"*_{layer_name}{layer_type}"
            )

            # Test if data was returned
            if len(paths) == 0:
                raise Exception(
                    f"No rasters found for grid cell {study_area} "
                    f"(raster version '{raster_version}'). Verify that "
                    f"`raster.py` has been run "
                    "for this grid cell."
                )

            # Create variable used for time axis
            time_var = xr.Variable("year", [int(i.split("/")[-1][0:4]) for i in paths])

            # Import data
            layer_da = xr.concat([xr.open_rasterio(i) for i in paths], dim=time_var)
            layer_da.name = f"{layer_name}"

            # Append to file
            da_list.append(layer_da)

        # Combine into a single dataset and restrict to start and end year
        layer_ds = xr.merge(da_list).squeeze("band", drop=True)
        layer_ds = layer_ds.sel(year=slice(str(start_year), str(end_year)))

        # Append to list
        ds_list.append(layer_ds)

    return ds_list


def ocean_masking(ds, ocean_da, connectivity=1, dilation=None):
    """
    Identifies ocean by selecting regions of water that overlap
    with ocean pixels. This region can be optionally dilated to
    ensure that the sub-pixel algorithm has pixels on either side
    of the water index threshold.

    Parameters:
    -----------
    ds : xarray.DataArray
        An array containing True for land pixels, and False for water.
        This can be obtained by thresholding a water index
        array (e.g. MNDWI < 0).
    ocean_da : xarray.DataArray
        A supplementary static dataset used to separate ocean waters
        from other inland water. The array should contain values of 1
        for high certainty ocean pixels, and 0 for all other pixels
        (land, inland water etc). For Australia, we use the  Geodata
        100K coastline dataset, rasterized as the "geodata_coast_100k"
        product on the DEA datacube.
    connectivity : integer, optional
        An integer passed to the 'connectivity' parameter of the
        `skimage.measure.label` function.
    dilation : integer, optional
        The number of pixels to dilate ocean pixels to ensure than
        adequate land pixels are included for subpixel waterline
        extraction. Defaults to None.

    Returns:
    --------
    ocean_mask : xarray.DataArray
        An array containing the a mask consisting of identified ocean
        pixels as True.
    """

    # Update `ocean_da` to mask out any pixels that are land in `ds` too
    ocean_da = ocean_da & (ds != 1)

    # First, break all time array into unique, discrete regions/blobs.
    # Fill NaN with 1 so it is treated as a background pixel
    blobs = xr.apply_ufunc(label, ds.fillna(1), 1, False, connectivity)

    # For each unique region/blob, use region properties to determine
    # whether it overlaps with a water feature from `water_mask`. If
    # it does, then it is considered to be directly connected with the
    # ocean; if not, then it is an inland waterbody.
    ocean_mask = blobs.isin(
        [i.label for i in regionprops(blobs.values, ocean_da.values) if i.max_intensity]
    )

    # Dilate mask so that we include land pixels on the inland side
    # of each shoreline to ensure contour extraction accurately
    # seperates land and water spectra
    if dilation:
        ocean_mask = xr.apply_ufunc(binary_dilation, ocean_mask, disk(dilation))

    return ocean_mask


def coastal_masking(ds, ocean_da, buffer=33):
    """
    Creates a symmetrical buffer around the land-water boundary
    in a input boolean array. This is used to create a study area
    mask that is focused on the coastal zone, excluding inland or
    deeper ocean pixels.

    Parameters:
    -----------
    ds : xarray.DataArray
        A single time-step boolean array containing True for land
        pixels, and False for water.
    ocean_da : xarray.DataArray
        A supplementary static dataset used to separate ocean waters
        from other inland water. The array should contain values of 1
        for high certainty ocean pixels, and 0 for all other pixels
        (land, inland water etc). For Australia, we use the  Geodata
        100K coastline dataset, rasterized as the "geodata_coast_100k"
        product on the DEA datacube.
    buffer : integer, optional
        The number of pixels to buffer the land-water boundary in
        each direction.

    Returns:
    --------
    coastal_mask : xarray.DataArray
        An array containing True within `buffer_pixels` of the
        land-water boundary, and False everywhere else.
    """

    def _coastal_buffer(ds, buffer):
        """Generate coastal buffer from ocean-land boundary"""
        buffer_ocean = binary_dilation(ds, buffer)
        buffer_land = binary_dilation(~ds, buffer)
        return buffer_ocean & buffer_land

    # Identify ocean pixels based on overlap with the Geodata
    # 100K coastline dataset
    all_time_ocean = ocean_masking(ds, ocean_da)

    # Generate coastal buffer from ocean-land boundary
    coastal_mask = xr.apply_ufunc(
        _coastal_buffer, all_time_ocean, disk(buffer), dask="parallelized"
    )

    # Return coastal mask as 1, and land pixels as 2
    return coastal_mask.where(coastal_mask, ~all_time_ocean * 2)


def temporal_masking(ds):
    """
    Create a temporal mask by identifying land pixels with a direct
    spatial connection (e.g. contiguous) to land pixels in either the
    previous or subsequent timestep.

    This is used to clean up noisy land pixels (e.g. caused by clouds,
    white water, sensor issues), as these pixels typically occur
    randomly with no relationship to the distribution of land in
    neighbouring timesteps. True land, however, is likely to appear
    in proximity to land before or after the specific timestep.

    Parameters:
    -----------
    ds : xarray.DataArray
        A multi-temporal array containing True for land pixels, and
        False for water.

    Returns:
    --------
    temporal_mask : xarray.DataArray
        A multi-temporal array array containing True for pixels
        located within the `dilation` distance of land in at least
        one neighbouring timestep.
    """

    def _noncontiguous(labels, intensity):

        # For each blob of land, obtain whether it intersected with land in
        # any neighbouring timestep
        region_props = regionprops(labels.values, intensity_image=intensity.values)
        contiguous = [i.label for i in region_props if i.max_intensity == 0]

        # Filter array to only contiguous land
        noncontiguous_array = np.isin(labels, contiguous)

        # Return as xr.DataArray
        return xr.DataArray(
            ~noncontiguous_array, coords=labels.coords, dims=labels.dims
        )

    # Label independent groups of pixels in each timestep in the array
    labelled_ds = xr.apply_ufunc(label, ds, None, 0, dask="parallelized").rename(
        "labels"
    )

    # Check if a pixel was neighboured by land in either the
    # previous or subsequent timestep by shifting array in both directions
    masked_neighbours = (
        (ds.shift(year=-1, fill_value=False) | ds.shift(year=1, fill_value=False))
        .astype(int)
        .rename("neighbours")
    )

    # Merge both into an xr.Dataset
    label_neighbour_ds = xr.merge([labelled_ds, masked_neighbours])

    # Apply continguity test to each year to obtain pixels that are
    # contiguous (spatially connected to) to land in the previous or subsequent timestep
    temporal_mask = label_neighbour_ds.groupby("year").apply(
        lambda x: _noncontiguous(labels=x.labels, intensity=x.neighbours)
    )

    return temporal_mask


def _create_mask(raster_mask, sieve_size, crs):
    """
    Clean and dilate an annual raster produced by `certainty_masking`,
    then vectorize into a dictionary of vector features that are
    taken as an input by `contour_certainty`.
    """

    # Clean mask by sieving to merge small areas of pixels into
    # their neighbours.
    sieved = xr.apply_ufunc(sieve, raster_mask, sieve_size)

    # Apply greyscale dilation to expand masked pixels and
    # err on the side of overclassifying certainty issues
    dilated = xr.apply_ufunc(dilation, sieved, disk(3))

    # Vectorise
    vector_mask = xr_vectorize(
        dilated,
        crs=crs,
        attribute_col="certainty",
    )

    # Dissolve column, fix geometry and rename classes
    vector_mask = vector_mask.dissolve("certainty")
    vector_mask["geometry"] = vector_mask.geometry.buffer(0)
    vector_mask = vector_mask.rename(
        {0: "good", 1: "unstable data", 2: "insufficient data"}
    )

    return (str(raster_mask.year.item()), vector_mask)


def certainty_masking(yearly_ds, obs_threshold=5, stdev_threshold=0.3, sieve_size=128):
    """
    Generate annual vector polygon masks containing information
    about the certainty of each extracted shoreline feature.
    These masks are used to assign each shoreline feature with
    important certainty information to flag potential issues with
    the data.

    Parameters:
    -----------
    yearly_ds : xarray.Dataset
        An `xarray.Dataset` containing annual DEA Coastlines
        rasters.
    obs_threshold : int, optional
        The minimum number of post-gapfilling Landsat observations
        required for an extracted shoreline to be considered good
        quality. Annual shorelines based on low numbers of
        observations can be noisy due to the influence of
        environmental noise like unmasked cloud, sea spray, white
        water etc. Defaults to 5.
    stdev_threshold : float, optional
        The maximum MNDWI standard deviation required for a
        post-gapfilled Landsat observation to be considered good
        quality. Annual shorelines based on MNDWI with a high
        standard deviation represent unstable data, which can
        indicate that the tidal modelling process did not adequately
        remove the influence of tide. For more information,
        refer to Bishop-Taylor et al. 2021
        (https://doi.org/10.1016/j.rse.2021.112734).
        Defaults to 0.3.
    sieve_size : int, optional
        To reduce the complexity of the output masks, they are
        first cleaned using `rasterio.features.sieve` to replace
        small areas of pixels with the values of their larger
        neighbours. This parameter sets the minimum polygon size
        to retain in this process. Defaults to 128.

    Returns:
    --------
    vector_masks : dictionary of geopandas.GeoDataFrames
        A dictionary with year (as an str) as the key, and vector
        data as a `geopandas.GeoDataFrame` for each year in the
        analysis.
    """

    from concurrent.futures import ProcessPoolExecutor
    from itertools import repeat

    # Identify problematic pixels
    high_stdev = yearly_ds["stdev"] > stdev_threshold
    low_obs = yearly_ds["count"] < obs_threshold

    # Create raster mask with values of 0 for good data, values of
    # 1 for unstable data, and values of 2 for insufficient data.
    raster_mask = high_stdev.where(~low_obs, 2).astype(np.int16)

    # Process in parallel
    with ProcessPoolExecutor() as executor:

        # Apply func in parallel, repeating params for each iteration
        groups = [group for (i, group) in raster_mask.groupby("year")]
        to_iterate = (
            groups,
            *(repeat(i, len(groups)) for i in [sieve_size, yearly_ds.odc.crs]),
        )
        vector_masks = dict(executor.map(_create_mask, *to_iterate), total=len(groups))

    return vector_masks


def contour_certainty(contours_gdf, certainty_masks):
    """
    Assigns a new certainty column to each annual shoreline feature
    to identify features affected by:

    1) Low satellite observations: annual shorelines based on less than
       5 annual observations after gapfilling
    2) Unstable MNDWI composites (potentially indicating tidal modelling
       issues): annual shorelines with MNDWI standard deviation > 0.3

    Parameters:
    -----------
    contours_gdf : geopandas.GeoDataFrame
        A `geopandas.GeoDataFrame` containing annual shorelines.
    certainty_masks : dictionary
        A dictionary of annual certainty mask vector features, as
        generated by `coastlines.vector.contours_preprocess`.

    Returns:
    --------
    contours_gdf : geopandas.GeoDataFrame
        A `geopandas.GeoDataFrame` containing annual shorelines with
        a new "certainty" column/field.
    """

    # Loop through each annual shoreline and attribute data with certainty
    out_list = []
    for year, _ in contours_gdf.iterrows():

        # Extract year
        contour_gdf = contours_gdf.loc[[year]]

        # Assign each shoreline segment with attributes from certainty mask
        contour_gdf = contour_gdf.overlay(
            certainty_masks[year].reset_index(), how="intersection"
        )

        # Set year field and use as index
        contour_gdf["year"] = year
        contour_gdf = contour_gdf.set_index("year")
        out_list.append(contour_gdf)

    # Combine into a single dataframe
    contours_gdf = pd.concat(out_list).sort_index()

    # Finally, set all 1991 and 1992 coastlines north of -23 degrees
    # latitude to 'uncertain' due to Mt Pinatubo aerosol issue
    pinatubo_lat = (contours_gdf.centroid.to_crs("EPSG:4326").y > -23) & (
        contours_gdf.index.isin(["1991", "1992"])
    )
    contours_gdf.loc[pinatubo_lat, "certainty"] = "aerosol issues"

    return contours_gdf


def contours_preprocess(
    yearly_ds,
    gapfill_ds,
    water_index,
    index_threshold,
    buffer_pixels=33,
    mask_temporal=True,
    mask_modifications=None,
    debug=False,
):
    """
    Prepares and preprocesses DEA Coastlines raster data to
    restrict the analysis to coastal shorelines, and extract data
    that is used to assess the certainty of extracted shorelines.

    This function:

    1) Identifies areas affected by either unstable composites or low data
    2) Fills low data areas in annual layers with three-year gapfill
    3) Computes an all-time coastal mask based on the observed timeseries
       of water and land pixels, after first optionally cleaning the data
       using land cover data, NDWI values and a temporal contiguity mask
    4) Identifies pixels directly connected to the ocean in each annual
       timestep, and uses this to remove inland waterbodies from the analyis

    Parameters:
    -----------
    yearly_ds : xarray.Dataset
        An `xarray.Dataset` containing annual DEA Coastlines rasters.
    gapfill_ds : xarray.Dataset
        An `xarray.Dataset` containing three-year gapfill DEA Coastlines
        rasters.
    water_index : string
        A string giving the name of the water index included in the
        annual and gapfill datasets (e.g. 'mndwi').
    index_threshold : float
        A float giving the water index threshold used to separate land
        and water (e.g. 0.00).
    buffer_pixels : int, optional
        The number of pixels by which to buffer the all time shoreline
        detected by this function to produce an overall coastal buffer.
        The default is 33 pixels, which at 30 m Landsat resolution
        produces a coastal buffer with a radius of approximately 1000 m.
    mask_temporal : bool, optional
        Whether to apply a temporal contiguity mask by identifying land
        pixels with a direct spatial connection (e.g. contiguous) to
        land pixels in either the previous or subsequent timestep. This is
        used to clean up noisy land pixels (e.g. caused by clouds,
        white water, sensor issues), as these pixels typically occur
        randomly with no relationship to the distribution of land in
        neighbouring timesteps. True land, however, is likely to appear
        in proximity to land before or after the specific timestep.
        Defaults to True.
    mask_modifications : geopandas.GeoDataFrame, optional
        An optional polygon dataset including features to remove or add
        to the all-time coastal mask. This should include a column/field
        named 'type' that contains two possible values:
            - 'add': features to add to the coastal mask (e.g. for
                     including areas of missing shorelines that were
                     previously removed by the coastal mask)
            - 'remove': features to remove from the coastal mask (e.g.
                        areas of non-coastal rivers or estuaries,
                        irrigated fields or aquaculture that you wish
                        to exclude from the analysis)
    debug : boolean, optional
        Whether to return all intermediate layers for troubleshooting.

    Returns:
    --------
    masked_ds : xarray.Dataset
        A dataset containing water index data for each annual timestep
        that has been masked to the coastal zone. This can then be used
        as an input to subpixel waterline extraction.
    certainty_masks : dict
        A dictionary containg one `geopandas.GeoDataFrame` for each year
        in the time period, with polygons identifying any potentially
        problematic region. This is used to assign each output shoreline
        with a certainty column.
    """

    # Remove low obs pixels and replace with 3-year gapfill
    combined_ds = yearly_ds.where(yearly_ds["count"] > 5, gapfill_ds)

    # Set any pixels with only one observation to NaN, as these are
    # extremely vulnerable to noise
    combined_ds = combined_ds.where(yearly_ds["count"] > 1)

    # Apply water index threshold and re-apply nodata values
    nodata = combined_ds[water_index].isnull()
    thresholded_ds = combined_ds[water_index] < index_threshold
    thresholded_ds = thresholded_ds.where(~nodata)

    # Compute temporal mask that restricts the analysis to land pixels
    # with a direct spatial connection (e.g. contiguous) to land pixels
    # in the previous or subsequent timestep. Set any pixels outside
    # mask to 0 to represent water
    temporal_mask = temporal_masking(thresholded_ds == 1)
    thresholded_ds = thresholded_ds.where(temporal_mask)

    # Create all time layer by identifying pixels that are land in at
    # least 15% of valid observations; this is used to generate an
    # all-of-time buffered coastal study area that constrains the Coastlines
    # analysis to the coastal zone
    all_time = thresholded_ds.mean(dim="year") > 0.5  # 0.15

    # Remove narrow river and stream features from all time layer using
    # the `black_tophat` transform. So narrow channels between islands
    # and the mainland aren't mistaken for rivers, first apply `sieve`
    # to any connected land pixels (i.e. islands) smaller than 5 km^2
    island_size = int(5000000 / (30 * 30))  # 5 km^2
    sieved = xr.apply_ufunc(sieve, all_time.astype(np.int16), island_size)
    rivers = xr.apply_ufunc(black_tophat, (sieved & all_time), disk(5)) == 1
    all_time_clean = all_time.where(~rivers, True)

    # Create a river mask by eroding the all time layer to clip out river
    # and estuary mouths, then expanding river features to include stream
    # banks and account for migrating rivers
    all_time_eroded = xr.apply_ufunc(binary_erosion, all_time_clean, disk(10))
    rivers = rivers.where(all_time_eroded, False)
    river_mask = ~xr.apply_ufunc(binary_dilation, rivers, disk(5))

    # Load Geodata 100K coastal layer to use to separate ocean waters from
    # other inland waters. This product has values of 0 for ocean waters,
    # and values of 1 and 2 for mainland/island pixels. We extract ocean
    # pixels (value 0), then erode these by 10 pixels to ensure we only
    # use high certainty deeper water ocean regions for identifying ocean
    # pixels in our satellite imagery
    dc = datacube.Datacube()
    geodata_da = dc.load(
        product="geodata_coast_100k",
        like=combined_ds.odc.geobox.compat,
    ).land.squeeze("time")
    ocean_da = xr.apply_ufunc(binary_erosion, geodata_da == 0, disk(10))

    # Use all time and Geodata 100K data to produce the buffered coastal
    # study area. Morphological closing helps to "close" the entrances of
    # estuaries and rivers, removing them from the analysis.
    coastal_mask = coastal_masking(
        ds=all_time_clean, ocean_da=ocean_da, buffer=buffer_pixels
    )

    # The output of `coastal_masking` contains values of 2 that represent
    # pixels inland of the coastal buffer. Use this to create
    # a mask of inland regions:
    inland_mask = coastal_mask != 2

    # Generate individual annual masks by selecting only water pixels that
    # are directly connected to the ocean in each yearly timestep
    annual_mask = (
        # Treat both 1s and NaN pixels as land (i.e. True)
        (thresholded_ds != 0)
        # Mask out rivers and inland regions
        .where(river_mask & inland_mask)
        # Keep pixels directly connected to ocean in each timestep
        .groupby("year").map(
            func=ocean_masking,
            ocean_da=ocean_da,
            connectivity=1,
            dilation=3,
        )
    )

    # Finally, apply temporal and annual masks to our surface water
    # index data, then clip to the all-of-time coastal study area
    masked_ds = combined_ds[water_index].where(
        temporal_mask & annual_mask & coastal_mask
    )

    # Generate annual vector polygon masks containing information
    # about the certainty of each shoreline feature
    certainty_masks = certainty_masking(combined_ds, stdev_threshold=0.3)

    # Return all intermediate layers if debug=True
    if debug:
        return (
            masked_ds,
            certainty_masks,
            all_time,
            all_time_clean,
            river_mask,
            ocean_da,
            coastal_mask,
            inland_mask,
            thresholded_ds,
            annual_mask,
        )

    else:
        return masked_ds, certainty_masks


def points_on_line(gdf, index, distance=30):
    """
    Generates evenly-spaced point features along a specific line feature
    in a `geopandas.GeoDataFrame`.

    Parameters:
    -----------
    gdf : geopandas.GeoDataFrame
        A `geopandas.GeoDataFrame` containing line features with an
        index and CRS.
    index : string or int
        An value giving the index of the line to generate points along
    distance : integer or float, optional
        A number giving the interval at which to generate points along
        the line feature. Defaults to 30, which will generate a point
        at every 30 metres along the line.

    Returns:
    --------
    points_gdf : geopandas.GeoDataFrame
        A `geopandas.GeoDataFrame` containing point features at every
        `distance` along the selected line.

    """

    # Select individual line to generate points along
    line_feature = gdf.loc[[index]].geometry

    # If multiple features are returned, take unary union
    if line_feature.shape[0] > 0:
        line_feature = line_feature.unary_union
    else:
        line_feature = line_feature.iloc[0]

    # Generate points along line and convert to geopandas.GeoDataFrame
    points_line = [
        line_feature.interpolate(i)
        for i in range(0, int(line_feature.length), distance)
    ]
    points_gdf = gpd.GeoDataFrame(geometry=points_line, crs=gdf.crs)

    return points_gdf


def annual_movements(
    points_gdf, contours_gdf, yearly_ds, baseline_year, water_index, max_valid_dist=1000
):
    """
    For each rate of change point along the baseline annual coastline,
    compute the distance to the nearest point on all neighbouring annual
    coastlines and add this data as new fields in the dataset.

    Distances are assigned a directionality (negative = located inland,
    positive = located sea-ward) by sampling water index values from the
    underlying DEA Coastlines rasters to determine if a coastline was
    located in wetter or drier terrain than the baseline coastline.

    Parameters:
    -----------
    points_gdf : geopandas.GeoDataFrame
        A `geopandas.GeoDataFrame` containing rates of change points.
    contours_gdf : geopandas.GeoDataFrame
        A `geopandas.GeoDataFrame` containing annual coastlines.
    yearly_ds : xarray.Dataset
        An `xarray.Dataset` containing annual DEA CoastLines rasters.
    baseline_year : string
        A string giving the year used as the baseline when generating
        the rates of change points dataset. This is used to load DEA
        CoastLines water index rasters to calculate change
        directionality.
    water_index : string
        A string giving the water index used in the analysis. This is
        used to load DEA CoastLines water index rasters to calculate
        change directionality.
    max_valid_dist : int or float, optional
        Any annual distance greater than this distance will be set
        to `np.nan`.

    Returns:
    --------
    points_gdf : geopandas.GeoDataFrame
        A `geopandas.GeoDataFrame` containing rates of change points
        with added 'dist_*' attribute columns giving the distance to
        each annual coastline from the baseline. Negative values
        indicate that an annual coastline was located inland of the
        baseline; positive values indicate the coastline was located
        towards the ocean.
    """

    def _point_interp(points, array, **kwargs):
        points_gs = gpd.GeoSeries(points)
        x_vals = xr.DataArray(points_gs.x, dims="z")
        y_vals = xr.DataArray(points_gs.y, dims="z")
        return array.interp(x=x_vals, y=y_vals, **kwargs)

    # Get array of water index values for baseline time period
    baseline_array = yearly_ds[water_index].sel(year=int(baseline_year))

    # Copy baseline point geometry to new column in points dataset
    points_gdf["p_baseline"] = points_gdf.geometry

    # Years to analyse
    years = contours_gdf.index.unique().values

    # Iterate through all comparison years in contour gdf
    for comp_year in years:

        # Set comparison contour
        comp_contour = contours_gdf.loc[[comp_year]].geometry.iloc[0]

        # Find nearest point on comparison contour, and add these to points dataset
        points_gdf[f"p_{comp_year}"] = points_gdf.apply(
            lambda x: nearest_points(x.p_baseline, comp_contour)[1], axis=1
        )

        # Compute distance between baseline and comparison year points and add
        # this distance as a new field named by the current year being analysed
        distances = points_gdf.apply(
            lambda x: x.geometry.distance(x[f"p_{comp_year}"]), axis=1
        )

        # Set any value over X m to NaN, and drop any points with
        # less than 50% valid observations
        points_gdf[f"dist_{comp_year}"] = distances.where(distances < max_valid_dist)

        # Extract comparison array containing water index values for the
        # current year being analysed
        comp_array = yearly_ds[water_index].sel(year=int(comp_year))

        # Sample water index values for baseline and comparison points
        points_gdf["index_comp_p1"] = _point_interp(
            points_gdf["p_baseline"], comp_array
        )
        points_gdf["index_baseline_p2"] = _point_interp(
            points_gdf[f"p_{comp_year}"], baseline_array
        )

        # Compute change directionality (positive = located towards the
        # ocean; negative = located inland)
        points_gdf["loss_gain"] = np.where(
            points_gdf.index_baseline_p2 > points_gdf.index_comp_p1, 1, -1
        )

        # Ensure NaNs are correctly propagated (otherwise, X > NaN
        # will return False, resulting in an incorrect land-ward direction)
        is_nan = points_gdf[["index_comp_p1", "index_baseline_p2"]].isna().any(axis=1)
        points_gdf["loss_gain"] = points_gdf["loss_gain"].where(~is_nan)

        # Multiply distance to set change to negative, positive or NaN
        points_gdf[f"dist_{comp_year}"] = (
            points_gdf[f"dist_{comp_year}"] * points_gdf.loss_gain
        )

        # Calculate compass bearing from baseline to comparison point;
        # first we need our points in lat-lon
        lat_lon = points_gdf[["p_baseline", f"p_{comp_year}"]].apply(
            lambda x: gpd.GeoSeries(x, crs=points_gdf.crs).to_crs("EPSG:4326")
        )

        geodesic = pyproj.Geod(ellps="WGS84")
        bearings = geodesic.inv(
            lons1=lat_lon.iloc[:, 0].values.x,
            lats1=lat_lon.iloc[:, 0].values.y,
            lons2=lat_lon.iloc[:, 1].values.x,
            lats2=lat_lon.iloc[:, 1].values.y,
        )[0]

        # Add bearing as a new column after first restricting
        # angles between 0 and 180 as we are only interested in
        # the overall axis of our points e.g. north-south
        points_gdf[f"bearings_{comp_year}"] = bearings % 180

    # Calculate mean and standard deviation of angles
    points_gdf["angle_mean"] = (
        points_gdf.loc[:, points_gdf.columns.str.contains("bearings_")]
        .apply(lambda x: circmean(x, high=180), axis=1)
        .round(0)
        .astype(int)
    )
    points_gdf["angle_std"] = (
        points_gdf.loc[:, points_gdf.columns.str.contains("bearings_")]
        .apply(lambda x: circstd(x, high=180), axis=1)
        .round(0)
        .astype(int)
    )

    # Keep only required columns
    to_keep = points_gdf.columns.str.contains("dist|geometry|angle")
    points_gdf = points_gdf.loc[:, to_keep]
    points_gdf = points_gdf.assign(**{f"dist_{baseline_year}": 0.0})
    points_gdf = points_gdf.round(2)

    return points_gdf


def outlier_mad(points, thresh=3.5):
    """
    Use robust Median Absolute Deviation (MAD) outlier detection
    algorithm to detect outliers. Returns a boolean array with True if
    points are outliers and False otherwise.

    Parameters:
    -----------
    points :
        An n-observations by n-dimensions array of observations
    thresh :
        The modified z-score to use as a threshold. Observations with a
        modified z-score (based on the median absolute deviation) greater
        than this value will be classified as outliers.

    Returns:
    --------
    mask :
        A n-observations-length boolean array.

    References:
    ----------
    Source: https://github.com/joferkington/oost_paper_code/blob/master/utilities.py

    Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
    Handle Outliers", The ASQC Basic References in Quality Control:
    Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    """
    if len(points.shape) == 1:
        points = points[:, None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median) ** 2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh


def outlier_ransac(xy_df, **kwargs):
    """
    Use the RANSAC (RANdom SAmple Consensus) algorithm to
    robustly identify outliers. Returns a boolean array with True if
    points are outliers and False otherwise.

    Parameters:
    -----------
    points :
        An n-observations by n-dimensions array of observations
    **kwargs :
        Any parameters to pass to
        `sklearn.linear_model.RANSACRegressor`

    Returns:
    --------
    mask :
        A n-observations-length boolean array.
    """

    from sklearn import linear_model

    # X and y inputs
    X = xy_df[:, 0].reshape(-1, 1)
    y = xy_df[:, 1].reshape(-1, 1)

    # Robustly fit linear model with RANSAC algorithm
    ransac = linear_model.RANSACRegressor(**kwargs)
    ransac.fit(X, y)
    inlier_mask = ransac.inlier_mask_
    outlier_mask = np.logical_not(inlier_mask)

    return outlier_mask


def change_regress(
    y_vals,
    x_vals,
    x_labels,
    threshold=3.5,
    detrend_params=None,
    slope_var="slope",
    interc_var="intercept",
    pvalue_var="pvalue",
    stderr_var="stderr",
    outliers_var="outliers",
):
    """
    For a given row in a `pandas.DataFrame`, apply linear regression to
    data values (as y-values) and a corresponding sequence of x-values,
    and return 'slope', 'intercept', 'pvalue', and 'stderr' regression
    parameters.

    Before computing the regression, outliers are identified using a
    robust Median Absolute Deviation (MAD) outlier detection algorithm,
    and excluded from the regression. A list of these outliers will be
    recorded in the output 'outliers' variable.

    Parameters:
    -----------
    x_vals, y_vals : list of numeric values, or nd.array
        A sequence of values to use as the x and y variables
    x_labels : list
        A sequence of strings corresponding to each value in `x_vals`.
        This is used to label any observations that are flagged as
        outliers (often, this can simply be set to the same list
        provided to `x_vals`).
    threshold : float, optional
        The modified z-score to use as a threshold for detecting
        outliers using the MAD algorithm. Observations with a modified
        z-score (based on the median absolute deviation) greater
        than this value will be classified as outliers.
    detrend_params : optional
        Not currently used
    slope, interc_var, pvalue_var, stderr_var : strings, optional
        Strings giving the names to use for each of the output
        regression variables.
    outliers_var : string, optional
        String giving the name to use for the output outlier variable.

    Returns:
    --------
    mask :
        A `pandas.Series` containing regression parameters and lists
        of outliers.

    """

    # Drop invalid NaN rows
    xy_df = np.vstack([x_vals, y_vals]).T
    valid_bool = ~np.isnan(xy_df).any(axis=1)
    xy_df = xy_df[valid_bool]
    valid_labels = x_labels[valid_bool]

    # If detrending parameters are provided, apply these to the data to
    # remove the trend prior to running the regression
    if detrend_params:
        xy_df[:, 1] = xy_df[:, 1] - (
            detrend_params[0] * xy_df[:, 0] + detrend_params[1]
        )

    # Remove outliers using MAD
    outlier_bool = outlier_mad(xy_df, thresh=threshold)
    # outlier_bool = outlier_ransac(xy_df)
    xy_df = xy_df[~outlier_bool]
    valid_labels = valid_labels[~outlier_bool]

    # Create string of all outliers and invalid NaN rows
    outlier_set = set(x_labels) - set(valid_labels)
    outlier_str = " ".join(map(str, sorted(outlier_set)))

    # Compute linear regression
    lin_reg = linregress(x=xy_df[:, 0], y=xy_df[:, 1])

    # Return slope, p-values and list of outlier years excluded from regression
    results_dict = {
        slope_var: np.round(lin_reg.slope, 3),
        interc_var: np.round(lin_reg.intercept, 3),
        pvalue_var: np.round(lin_reg.pvalue, 3),
        stderr_var: np.round(lin_reg.stderr, 3),
        outliers_var: outlier_str,
    }

    return pd.Series(results_dict)


def calculate_regressions(points_gdf):
    """
    For each rate of change point along the baseline annual coastline,
    compute linear regression rates of change against both time and
    climate indices.

    Regressions are computed after removing outliers to ensure robust
    results.

    Parameters:
    -----------
    points_gdf : geopandas.GeoDataFrame
        A `geopandas.GeoDataFrame` containing rates of change points
        with 'dist_*' annual movement/distance data.

    Returns:
    --------
    points_gdf : geopandas.GeoDataFrame
        A `geopandas.GeoDataFrame` containing rates of change points
        with additional attribute columns:

            'rate_*':  Slope of the regression
            'sig_*':   Significance of the regression
            'se_*':    Standard error of the  regression
            'outl_*':  A list of any outlier years excluded from the
                       regression
    """

    # Restrict data to years in datasets
    dist_years = points_gdf.columns[points_gdf.columns.str.contains("dist_")]
    x_years = dist_years.str.replace("dist_", "").astype(int)
    points_subset = points_gdf[dist_years]

    # Compute coastal change rates by linearly regressing annual
    # movements vs. time
    rate_out = points_subset.apply(
        lambda row: change_regress(
            y_vals=row.values.astype(float), x_vals=x_years, x_labels=x_years
        ),
        axis=1,
    )
    points_gdf[
        ["rate_time", "incpt_time", "sig_time", "se_time", "outl_time"]
    ] = rate_out

    # Copy slope and intercept into points_subset so they can be
    # used to temporally de-trend annual distances
    points_subset[["slope", "intercept"]] = rate_out[["slope", "intercept"]]

    # Custom sorting
    reg_cols = ["rate_time", "sig_time", "se_time", "outl_time"]

    return points_gdf.loc[
        :, [*reg_cols, *dist_years, "angle_mean", "angle_std", "geometry"]
    ]


def all_time_stats(x, col="dist_", initial_year=1988):
    """
    Apply any statistics that apply to the entire set of annual
    distance/movement values. This currently includes:

        valid_obs, valid_span : The number of valid (non-outlier)
             obervations, and the length of time in years between
             the first and last valid observation.
        sce: Shoreline Change Envelope (SCE). A measure of the maximum
             change or variability across all annual coastlines,
             calculated by computing the maximum distance between any
             two annual coastlines (excluding outliers).
        nsm: Net Shoreline Movement (NSM). The distance between the
             oldest and most recent annual shorelines (excluding
             outliers). Negative values indicate the shoreline retreated
             between the oldest and most recent shoreline; positive
             values indicate growth.
        max_year, min_year: The year that annual shorelines were at
             their maximum (i.e. located furthest towards the ocean) and
             their minimum (i.e. located furthest inland) respectively
             (excluding outliers).

    Parameters:
    -----------
    x : pandas.DataFrame row
        A single row of the annual rates of change `pandas.DataFrame`
        containg columns of annual distances from the baseline.
    col : string, optional
        A string giving the prefix used for all annual distance/
        movement values. The default is 'dist_'.
    initial_year : int, optional
        An optional integer giving the first year of data to use when
        calculating statistics. This can be useful when data from early
        in the satellite timeseries is less reliable than more recent
        data, e.g. in regions with sparse Landsat 5 satellite coverage.

    Returns:
    --------
    A `pandas.Series` containing new all time statistics.
    """

    # Select date columns only
    year_cols = x.index.str.contains(col)
    subset = x.loc[year_cols].astype(float)

    # Restrict to requested initial year
    subset.index = subset.index.str.lstrip("dist_").astype(int)
    subset = subset.loc[initial_year:]

    # Identify outlier years to drop from calculation
    to_drop = [int(i) for i in x.outl_time.split(" ") if len(i) > 0]
    subset_nooutl = subset.drop(to_drop, errors="ignore")

    # Calculate SCE range, NSM and max/min year
    # Since NSM is the most recent shoreline minus the oldest shoreline,
    # we can calculate this by simply inverting the 1988 distance value
    # (i.e. 0 - X) if it exists in the data
    stats_dict = {
        "valid_obs": subset_nooutl.shape[0],
        "valid_span": (subset_nooutl.index[-1] - subset_nooutl.index[0] + 1),
        "sce": subset_nooutl.max() - subset_nooutl.min(),
        "nsm": -(
            subset_nooutl.loc[initial_year] if initial_year in subset_nooutl else np.nan
        ),
        "max_year": subset_nooutl.idxmax(),
        "min_year": subset_nooutl.idxmin(),
    }

    return pd.Series(stats_dict)


def rocky_shoreline_flag(
    points_gdf,
    geomorphology_gdf,
    rocky_query="(Preds == 'Bedrock') and (Probs > 0.75)",
    max_distance=300,
):
    """
    Identifies rate of change points that are potentially being part
    of a rocky shoreline using geomorphology classification data.
    The input geomorphology data should contain attributes that can
    be queried to identify rocky shorelines.

    This can be important as deep shadows in areas of rocky shorelines
    can produce misleading impressions of coastal change (particularly
    in non-terrain corrected satellite data).

    Parameters:
    -----------
    points_gdf : geopandas.GeoDataFrame
        A `geopandas.GeoDataFrame` containing rates of change points.
    geomorphology_gdf : geopandas.GeoDataFrame
        A `geopandas.GeoDataFrame` containing geomorphology classification
        data that can be analysed to identify rocky shorelines.
    rocky_query : str, optional
        A string that can be passed as an input to the Pandas `.eval`
        method to filter features from `geomorphology_gdf` to identify
        rocky shorelines.
    max_distance : int, optional
        The maximum distance to join geomorphology features to each
        rate of change point. Defaults to 300 metres.

    Returns:
    --------
    pandas.Series
        A column with True if a point is likely to be a rocky shoreline;
        otherwise False.
    """

    # Classify geomorphology data to rocky shores or not using query
    geomorphology_gdf["rocky"] = geomorphology_gdf.eval(rocky_query)

    # Join classified geomorphology data to points if within max dist
    joined = gpd.sjoin_nearest(
        points_gdf,
        geomorphology_gdf[["rocky", "geometry"]],
        how="left",
        max_distance=300,
    )

    # Return boolean indicating whether point was rocky; take max of
    # each unique index value (i.e. True if there are both True and False)
    # to account for edge case where nearest geomorphology is the corner
    # of two vector features
    return (joined["rocky"] == True).groupby(joined.index).max()


def region_atttributes(gdf, region_gdf, attribute_col="TERRITORY1", rename_col=False):
    """
    Produces an attribute column for each rates of change point or
    annual shoreline in the dataset by spatially joining regions from an
    external vector file. This can be used, for example, to assign
    a "country" or "region" column to each spatial point or coastline.

    Parameters:
    -----------
    gdf : geopandas.GeoDataFrame
        A `geopandas.GeoDataFrame` containing rates of change points or
        annual shoreline vectors.
    region_gdf : geopandas.GeoDataFrame
        A `geopandas.GeoDataFrame` containing region data that will be
        spatially joined to feature in `gdf`.
    attribute_col : str or list, optional
        A string (or list of strings) providing the names of the
        attribute columns from `region_gdf` you wish to obtain.
    rename_col : str or list, optional
        An option string (or list of strings) giving new names for
        each attribute columns specified in `attribute_col`. If passing
        a list, ensure this is the same length as `attribute_col`.

    Returns:
    --------
    geopandas.GeoDataFrame
        The `geopandas.GeoDataFrame` provided to `gdf` with one or multiple
        additional attribute columns added from `region_gdf`.
    """

    # Wrap values in list if provided as strings
    if not isinstance(attribute_col, (list, tuple)):
        attribute_col = [attribute_col]
        rename_col = [rename_col]

    # Select subset of attribute columns to join (need to include
    # geometry as well to enable spatial join)
    region_subset = region_gdf[[*attribute_col, "geometry"]]

    # Rename columns if requested
    if rename_col:
        region_subset = region_subset.rename(
            dict(zip(attribute_col, rename_col)), axis=1
        )

    # Spatial join region data to points
    if gdf.iloc[0].geometry.type == "Point":
        joined_df = gdf.sjoin(region_subset, how="left").drop("index_right", axis=1)

    # Or if data is not points, use overlay (overlay removes index on
    # gdf1, so we need to reset to keep it as a columnm, then reapply)
    else:
        joined_df = gpd.overlay(
            gdf.reset_index(),
            region_subset,
            how="union",
            keep_geom_type=True,
        ).set_index(gdf.index.name)

    return joined_df


def vector_schema(gdf, default="float:8.2"):
    """
    Creates a schema passed to `gdf.to_file(schema=...)` that
    sets the dtype and precision of each DEA Coastlines column
    in `gdf`.

    This helps to greatly reduce the size of the output
    files on disk.

    Parameters:
    -----------
    gdf : geopandas.GeoDataFrame
        The input vector dataset containing one or more columns/
        attribute fields.
    default : string, optional
        An optional string giving the default dtype/precision of
        any column in `gdf` that does not exist in the list of
        custom fields below. Defaults to 'float:8.2' for a floating
        point column with two decimal places precision.

    Returns:
    --------
    schema_dict : dict
        A dictionary mapping all columns in `gdf` to a dtype/precision.
        This can be passed to `gdf.to_file(schema=...)`.
    """

    # Create default value for all fields
    schema_dict = {key: "float:8.2" for key in gdf.reset_index().columns}

    # Update default to custom values
    schema_dict.update(
        {
            # Rates of change
            "uid": "str:11",
            "sig_time": "float:8.3",
            "outl_time": "str:80",
            "angle_mean": "int:3",
            "angle_std": "int:3",
            "valid_obs": "int:4",
            "valid_span": "int:4",
            "max_year": "int:4",
            "min_year": "int:4",
            "certainty": "str:25",
            "id_primary": "str:10",
            # Annual shorelines only
            "year": "int:4",
            "tide_datum": "str:20",
            # Hotspots only
            "n": "int:6",
            "radius_m": "int:6",
            # WMS only
            "wms_conf": "float:8.1",
            "wms_grew": "int:1",
            "wms_retr": "int:1",
            "wms_sig": "int:1",
        }
    )

    # Filter to columns in data
    return {
        key: schema_dict[key] for key in gdf.reset_index().columns if key != "geometry"
    }


def generate_vectors(
    config,
    study_area,
    raster_version,
    vector_version,
    water_index,
    index_threshold,
    start_year,
    end_year,
    baseline_year,
    log=None,
):
    ###############################
    # Load DEA Coastlines rasters #
    ###############################

    if log is not None:
        log = configure_logging()

    log.info(f"Study area {study_area}: Starting vector generation")

    yearly_ds, gapfill_ds = load_rasters(
        path="data/interim/raster",
        raster_version=raster_version,
        study_area=study_area,
        water_index=water_index,
        start_year=start_year,
        end_year=end_year,
    )
    log.info(f"Study area {study_area}: Loaded rasters")

    # Create output vector folder using supplied vector version string;
    output_dir = f"data/interim/vector/{vector_version}/{study_area}_{vector_version}"
    os.makedirs(output_dir, exist_ok=True)

    ####################
    # Load vector data #
    ####################

    # Get bounding box to load data for
    bbox = gpd.GeoSeries(yearly_ds.odc.geobox.extent.geom, crs=yearly_ds.odc.crs)

    # Study area polygon
    gridcell_gdf = (
        gpd.read_file(config["Input files"]["grid_path"], bbox=bbox)
        .set_index("id")
        .to_crs(str(yearly_ds.odc.crs))
    )
    gridcell_gdf.index = gridcell_gdf.index.astype(int).astype(str)
    gridcell_gdf = gridcell_gdf.loc[[str(study_area)]]

    #     # Coastal mask modifications
    #     modifications_gdf = gpd.read_file(
    #         config["Input files"]["modifications_path"], bbox=bbox
    #     ).to_crs(str(yearly_ds.odc.crs))

    # Geomorphology dataset
    geomorphology_gdf = gpd.read_file(
        config["Input files"]["geomorphology_path"], bbox=bbox
    ).to_crs(str(yearly_ds.odc.crs))

    # Region attribute dataset
    region_gdf = gpd.read_file(
        config["Input files"]["region_attributes_path"], bbox=bbox
    ).to_crs(str(yearly_ds.odc.crs))

    ##############################
    # Extract shoreline contours #
    ##############################

    # Mask dataset to focus on coastal zone only
    masked_ds, certainty_masks = contours_preprocess(
        yearly_ds,
        gapfill_ds,
        water_index,
        index_threshold,
        buffer_pixels=33,
        mask_modifications=None,
    )

    # Extract annual shorelines
    contours_gdf = subpixel_contours(
        da=masked_ds,
        z_values=index_threshold,
        min_vertices=10,
        dim="year",
    ).set_index("year")

    if len(contours_gdf.index) == 0:
        raise ValueError(
            f"Study area {study_area}: Unable to extract any valid shorelines from raster data"
        )

    log.info(f"Study area {study_area}: Extracted shorelines from raster data")

    ######################
    # Compute statistics #
    ######################

    # Extract statistics modelling points along baseline shoreline
    try:

        points_gdf = points_on_line(contours_gdf, str(baseline_year), distance=30)
        log.info(f"Study area {study_area}: Extracted rates of change points")

    except KeyError:

        log.warning(
            f"Study area {study_area}: Baseline year {baseline_year} missing from annual shorelines; unable to extract rates of change points"
        )
        points_gdf = None

    # If any points exist in the dataset
    if points_gdf is not None and len(points_gdf) > 0:

        # Calculate annual coastline movements and residual tide heights
        # for every contour compared to the baseline year
        points_gdf = annual_movements(
            points_gdf,
            contours_gdf,
            yearly_ds,
            str(baseline_year),
            water_index,
            max_valid_dist=1200,
        )

        # Reindex to add any missing annual columns to the dataset
        points_gdf = points_gdf.reindex(
            columns=[
                "geometry",
                *[f"dist_{i}" for i in range(start_year, end_year + 1)],
                "angle_mean",
                "angle_std",
            ]
        )
        log.info(
            f"Study area {study_area}: Calculated distances to each annual shoreline"
        )

        # Calculate regressions
        points_gdf = calculate_regressions(points_gdf)
        log.info(f"Study area {study_area}: Calculated rates of change regressions")

        # Add count and span of valid obs, Shoreline Change Envelope
        # (SCE), Net Shoreline Movement (NSM) and Max/Min years
        stats_list = ["valid_obs", "valid_span", "sce", "nsm", "max_year", "min_year"]
        points_gdf[stats_list] = points_gdf.apply(
            lambda x: all_time_stats(x, initial_year=start_year), axis=1
        )
        log.info(f"Study area {study_area}: Calculated all of time statistics")

        # Add certainty column to flag points with:
        # - Baseline outlier: The baseline shoreline is itself flagged as an
        #   outlier, potentially resulting in inaccurate rates of change.
        # - Likely rocky shorelines: Rates of change can be unreliable in areas
        #   with steep rocky/bedrock shorelines due to terrain shadow.
        # - Extreme rate of change value (> 50 m per year change): these are more
        #   likely to reflect modelling issues than real-world coastal change
        # - High angular variability: the nearest shorelines for each year do not
        #   fall on an approximate line, making rates of change invalid
        # - Insufficient observations: less than 75% valid annual shorelines, which
        #   make the resulting rates of change more likely to be inaccurate
        rocky = [
            "Bedrock breakdown debris (cobbles/boulders)",
            "Boulder (rock) beach",
            "Cliff (>5m) (undiff)",
            "Colluvium (talus) undiff",
            "Flat boulder deposit (rock) undiff",
            "Hard bedrock shore",
            "Hard bedrock shore inferred",
            "Hard rock cliff (>5m)",
            "Hard rocky shore platform",
            "Rocky shore (undiff)",
            "Rocky shore platform (undiff)",
            "Sloping hard rock shore",
            "Sloping rocky shore (undiff)",
            "Soft `bedrock¿ cliff (>5m)",
            "Steep boulder talus",
            "Hard rocky shore platform",
        ]

        # Initialise certainty column with good values
        points_gdf["certainty"] = "good"

        # Flag points where the baseline shoreline is itself an outlier
        points_gdf.loc[
            points_gdf.outl_time.str.contains(str(baseline_year)), "certainty"
        ] = "baseline outlier"

        # Flag rocky shorelines
        points_gdf.loc[
            rocky_shoreline_flag(
                points_gdf,
                geomorphology_gdf,
                rocky_query=f"(INTERTD1_V in {rocky}) & (INTERTD2_V in {rocky + ['Unclassified']})",
            ),
            "certainty",
        ] = "likely rocky coastline"

        # Flag extreme rates of change
        points_gdf.loc[
            points_gdf.rate_time.abs() > 50, "certainty"
        ] = "extreme value (> 50 m)"

        # Flag points where change does not fall on a line
        points_gdf.loc[
            points_gdf.angle_std > 30, "certainty"
        ] = "high angular variability"

        # Flag shorelines with less than X valid shorelines
        valid_obs_thresh = int(points_gdf.columns.str.contains("dist_").sum() * 0.75)
        points_gdf.loc[
            points_gdf.valid_obs < valid_obs_thresh, "certainty"
        ] = "insufficient observations"

        log.info(f"Study area {study_area}: Calculated rate of change certainty flags")

        # Add region attributes
        points_gdf = region_atttributes(
            points_gdf, region_gdf, attribute_col="ID_Primary", rename_col="id_primary"
        )

        # Generate a geohash UID for each point and set as index
        uids = (
            points_gdf.geometry.to_crs("EPSG:4326")
            .apply(lambda x: gh.encode(x.y, x.x, precision=10))
            .rename("uid")
        )
        points_gdf = points_gdf.set_index(uids)

        log.info(f"Study area {study_area}: Added region attributes and geohash UIDs")

        ################
        # Export stats #
        ################

        #         if points_gdf is not None and len(points_gdf) > 0:

        # Clip stats to study area extent
        points_gdf = points_gdf[points_gdf.intersects(gridcell_gdf.geometry.item())]

        # Set output path
        stats_path = (
            f"{output_dir}/ratesofchange_{study_area}_"
            f"{vector_version}_{water_index}_{index_threshold:.2f}"
        )

        try:

            # Export to GeoJSON
            points_gdf.to_crs("EPSG:4326").to_file(
                f"{stats_path}.geojson",
                driver="GeoJSON",
            )

            # Export as ESRI shapefiles
            points_gdf.to_file(
                f"{stats_path}.shp",
                schema={
                    "properties": vector_schema(points_gdf),
                    "geometry": "Point",
                },
            )

        except ValueError:
            raise ValueError(
                f"Study area {study_area}: No vector points data to export after clipping to study area extent"
            )

    else:
        log.warning(f"Study area {study_area}: No rates of change points to process")

    #####################
    # Export shorelines #
    #####################

    # Assign certainty to shorelines based on underlying masks
    contours_gdf = contour_certainty(contours_gdf, certainty_masks)

    # Add tide datum details (this supports future addition of extra tide datums)
    contours_gdf["tide_datum"] = "0 m AMSL"

    # Add region attributes
    contours_gdf = region_atttributes(
        contours_gdf, region_gdf, attribute_col="ID_Primary", rename_col="id_primary"
    )

    # Set output path
    contour_path = (
        f"{output_dir}/annualshorelines_{study_area}_{vector_version}_"
        f"{water_index}_{index_threshold:.2f}"
    )

    # Clip annual shoreline contours to study area extent
    contours_gdf["geometry"] = contours_gdf.intersection(gridcell_gdf.geometry.item())

    try:

        # Export to GeoJSON
        contours_gdf.to_crs("EPSG:4326").to_file(
            f"{contour_path}.geojson", driver="GeoJSON"
        )

        # Export stats and contours as ESRI shapefiles
        contours_gdf.to_file(
            f"{contour_path}.shp",
            schema={
                "properties": vector_schema(contours_gdf),
                "geometry": ["MultiLineString", "LineString"],
            },
        )

    except ValueError:
        raise ValueError(
            f"Study area {study_area}: No vector shorelines data to export after clipping to study area extent"
        )

    log.info(f"Study area {study_area}: Output vector files written to {output_dir}")


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
    "gridcell that was previously used to generate raster "
    "files. This is used to identify the raster files that "
    "will be used as inputs for shoreline extraction, and "
    'should match a row in the "id" column of the provided '
    "analysis gridcell vector file.",
)
@click.option(
    "--raster_version",
    type=str,
    required=True,
    help="A unique string providing a name that was used "
    "to generate raster files. This is used to identify the "
    "raster files that will be used as inputs for shoreline "
    "extraction.",
)
@click.option(
    "--vector_version",
    type=str,
    help="A unique string proving a name that will be used "
    "for output vector directories and files. This allows "
    "multiple versions of vector files to be generated "
    "from the same input raster data, e.g. for testing "
    "different water index thresholds or indices. If "
    "not provided, this will default to the same string "
    'supplied to "--raster_version".',
)
@click.option(
    "--water_index",
    type=str,
    default="mndwi",
    help="A string giving the name of the computed water "
    "index to use for shoreline extraction. "
    'Defaults to "mndwi".',
)
@click.option(
    "--index_threshold",
    type=float,
    default=0.00,
    help="The water index threshold used to extract "
    "subpixel precision shorelines. Defaults to 0.00.",
)
@click.option(
    "--start_year",
    type=int,
    default=1988,
    help="The first annual shoreline to extract from the input raster data.",
)
@click.option(
    "--end_year",
    type=int,
    default=2021,
    help="The final annual shoreline to extract from the input raster data.",
)
@click.option(
    "--baseline_year",
    type=int,
    default=2021,
    help="The annual shoreline used as a baseline from "
    "which to generate the rates of change point statistics. "
    "This is typically the most recent annual shoreline in "
    "the dataset (i.e. the same as `--end_year`).",
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
def generate_vectors_cli(
    config_path,
    study_area,
    raster_version,
    vector_version,
    water_index,
    index_threshold,
    start_year,
    end_year,
    baseline_year,
    aws_unsigned,
    overwrite,
):

    log = configure_logging(f"Coastlines vector generation for study area {study_area}")

    # If no vector version is provided, copy raster version
    if vector_version is None:
        vector_version = raster_version

    # Test if study area has already been run by checking if run status file exists
    run_status_file = f"data/interim/vector/{vector_version}/{study_area}_{vector_version}/run_completed"
    output_exists = os.path.exists(run_status_file)

    # Skip if outputs exist but overwrite is False
    if output_exists and not overwrite:
        log.info(
            f"Study area {study_area}: Data exists but overwrite set to False; skipping."
        )
        sys.exit(0)

    # Load analysis params from config file
    config = load_config(config_path=config_path)

    # Do an opinionated configuration of S3
    configure_s3_access(cloud_defaults=True, aws_unsigned=aws_unsigned)

    # Run the code to generate vectors
    try:
        generate_vectors(
            config,
            study_area,
            raster_version,
            vector_version,
            water_index,
            index_threshold,
            start_year,
            end_year,
            baseline_year,
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
    generate_vectors_cli()
