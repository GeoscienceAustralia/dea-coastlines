#!/usr/bin/env python
# coding: utf-8

# This code combines individual datasets into continental DE Africa
# Coastlines layers:
#
#     * Combines output shorelines and rates of change statistics point
#       vectors into single continental datasets
#     * Aggregates this data to produce moving window hotspot datasets
#       that summarise coastal change at regional and continental scale.

import os
import sys
import click
import fiona
import geopandas as gpd
import pandas as pd
from rtree import index
from tqdm.auto import tqdm

# Import DEA Coastlines code
from dea_coastlines import vector


def points_in_poly(points, polygons):
    """
    Builds an optimised spatial index using `rtree` to determine the IDs
    of points that fall within each polygon ID.

    Parameters:
    -----------
    points :
        A `geopandas.GeoSeries` or iterable of `shapely` geometry points
    polygons :
        A `geopandas.GeoSeries` or iterable of `shapely` geometry polygons

    Returns:
    --------
    A dictionary identifying what point IDs fall within each polygon ID.
    the polygon.
    """

    # Create the R-tree index and store the features in it (bounding box)
    idx = index.Index()
    for pos, poly in enumerate(tqdm(polygons, desc='Building index')):
        idx.insert(pos, poly.bounds)

    # Iterate through points
    out_dict = {}
    for i, point in enumerate(tqdm(points, desc='Processing points')):
        poly_ids = [
            j for j in idx.intersection((point.coords[0]))
            if point.within(polygons[j])
        ]
        out_dict[i] = poly_ids

    return out_dict


def get_matching_data(key, stats_gdf, poly_points_dict, min_n=100):
    """
    Computes statistics based on all spatial points that intersect
    with a specific polygon. This is used to calculate moving
    window statistics after first creating a dictionary of polygons
    and intersecting points using `points_in_poly`.

    Parameters:
    -----------
    key :
        The ID of a polygon that will be used to select all
        intersecting points.
    stats_gdf : geopandas.GeoDataFrame
        A `geopandas.GeoDataFrame` containing multiple points which
        will be selected based on whether they match the points linked
        to the specific polygon in `poly_points_dict` requested by `key`.
    poly_points_dict : dict
        A dictionary identifying what point IDs fall within each polygon
        ID. An individual polygon is selected using `key`.
    min_n : int, optional
        If less than `min_n` points are returned for a given polygon,
        return None.

    Returns:
    --------
    A `pandas.Series` containing statistics based on the points that
    matched the polygon ID requested by `key`.
    """

    matching_points = stats_gdf.iloc[poly_points_dict[key]]

    if len(matching_points.index) > min_n:

        return pd.Series([
            matching_points.rate_time.mean(),
            len(matching_points.index)
        ])

    else:
        return pd.Series([None, None])


@click.command()
@click.option('--vector_version',
              type=str,
              required=True,
              help='A unique string proving a name that was used '
              'for output vector directories and files. This is used '
              'to identify the tiled annual shoreline and rates of '
              'change layers that will be combined into continental-'
              'scale layers.')
@click.option('--continental_version',
              type=str,
              help='A unique string proving a name that will be used '
              'for output continental-scale layers. This allows '
              'multiple versions of continental-scale layers to be '
              'generated from the same input vector data, e.g. for '
              'testing different hotspot of coastal change summary '
              'layers. If not provided, this will default to the '
              'string provided to "--vector_version".')
@click.option('--water_index',
              type=str,
              default='mndwi',
              help='The water index used to extract annual shorelines. '
              'Used to identify tiled annual shoreline and rates of '
              'change layers to combine into individual continental-'
              'scale layers. Defaults to "mndwi".')
@click.option('--index_threshold',
              default='0.00',
              type=str,
              help='The water index threshold used to extract annual '
              'shorelines. Used to identify tiled annual shoreline and '
              'rates of change layers to combine into individual '
              'continental-scale layers. Defaults to "0.00".')
@click.option('--shorelines',
              type=bool,
              default=True,
              help='A boolean indicating whether to combine tiled '
              'annual shorelines layers into a single continental-'
              'scale annual shorelines layer.')
@click.option('--ratesofchange',
              type=bool,
              default=True,
              help='A boolean indicating whether to combine tiled '
              'rates of change statistics layers into a single '
              'continental-scale rates of change statistics layer.')
@click.option('--hotspots',
              help='The distance (in metres) used to generate a '
              'coastal change hotspots summary layer. This controls '
              'the spacing of each summary point, and the radius used '
              'to aggregate rates of change statistics around each '
              'point. The default is `False` which will not generate a '
              'hotspots layer. If `True` is supplied instead of a '
              'distance value, hotspots will be generated using a '
              'default distance of 10000 m.')
@click.option('--baseline_year',
              type=str,
              default='2020',
              help='The annual shoreline used to generate the hotspot '
              'summary points. This is typically the most recent '
              'annual shoreline in the dataset.')
def continental_layers(vector_version, continental_version, water_index,
                       index_threshold, shorelines, ratesofchange, hotspots,
                       baseline_year):

    #################
    # Merge vectors #
    #################

    # If no continental version is provided, copy this from vector
    # version
    if continental_version is None:
        continental_version = vector_version
    output_dir = f'data/processed/{continental_version}'    
    os.makedirs(output_dir, exist_ok=True)

    # Setup input and output file paths
    shoreline_paths = f'data/interim/vector/{vector_version}/*/' \
                      f'annualshorelines_*_{vector_version}_' \
                      f'{water_index}_{index_threshold}.shp'
    ratesofchange_paths = f'data/interim/vector/{vector_version}/*/' \
                          f'ratesofchange_*_{vector_version}_' \
                          f'{water_index}_{index_threshold}.shp'
    continental_shorelines_path = f'{output_dir}/DEAfricaCoastlines_' \
                                  f'annualshorelines_{continental_version}.shp'
    continental_rates_path = f'{output_dir}/DEAfricaCoastlines_' \
                             f'ratesofchange_{continental_version}.shp'

    # Combine annual shorelines into a single continental layer
    if shorelines:
        print('Combining annual shorelines...')
        os.system(f'ogrmerge.py -o '
                  f'{continental_shorelines_path} {shoreline_paths} '
                  f'-single -overwrite_ds -t_srs EPSG:3577')

    # Combine rates of change stats points into single continental layer
    if ratesofchange:
        print('Combining rates of change statistics...')
        os.system(f'ogrmerge.py '
                  f'-o {continental_rates_path} {ratesofchange_paths} '
                  f'-single -overwrite_ds -t_srs EPSG:3577')

    # Generate hotspot points that provide regional/continental summary
    # of hotspots of coastal erosion and growth
    if hotspots:

        ###############################
        # Load DEA CoastLines vectors #
        ###############################

        print('Generating hotspots...')

        # If hotspots is provided as True rather than a distance, use
        # a default distance radius of 10000 m
        if hotspots is True:
            print('Using a default hotspot distance of 10000 m')
            hotspots = 10000
        hotspots = int(hotspots)

        # Load continental shoreline and rates of change data
        try:
            ratesofchange_gdf = gpd.read_file(continental_rates_path)
            shorelines_gdf = gpd.read_file(continental_shorelines_path)
        except fiona.errors.DriverError:
            raise FileNotFoundError(
                'Continental-scale annual shoreline and rates of '
                'change layers are required for hotspot generation. '
                'Try re-running this analysis with the following '
                'settings: `--shorelines True --ratesofchange True`.')

        # Set year index on coastlines
        shorelines_gdf = (shorelines_gdf.loc[
            shorelines_gdf.geometry.is_valid].set_index('year'))

        # Extract hotspot points
        hotspots_gdf = vector.points_on_line(shorelines_gdf,
                                             index=baseline_year,
                                             distance=hotspots)

        # Drop low observations (less than 10) from rates
        ratesofchange_gdf = ratesofchange_gdf.loc[
            ratesofchange_gdf.valid_obs > 10]
        ratesofchange_gdf = ratesofchange_gdf.reset_index(drop=True)

        # Set nonsignificant rates to 0 m / year
        ratesofchange_gdf.loc[ratesofchange_gdf.sig_time > 0.01,
                              'rate_time'] = 0

        # Clip to 50 m rates to remove extreme outliers
        ratesofchange_gdf['rate_time'] = ratesofchange_gdf.rate_time.clip(
            -50, 50)

        #####################
        # Generate hotspots #
        #####################

        # Generate dictionary of polygon IDs and corresponding points
        print('Identifying points in each polygon')
        poly_points_dict = points_in_poly(points=hotspots_gdf.geometry,
                                          polygons=ratesofchange_gdf.buffer(
                                              hotspots * 2))

        # Compute mean and number of obs for each polygon
        print('Calculating mean rates')
        hotspots_gdf[['rate_time', 'n']] = hotspots_gdf.apply(
            lambda row: get_matching_data(row.name, 
                                          ratesofchange_gdf,
                                          poly_points_dict,
                                          min_n=hotspots / 25), axis=1)

        # Export hotspots to file
        hotspots_gdf.to_file(f'{output_dir}/DEAfricaCoastlines_hotspots_'
                             f'{continental_version}_{hotspots}.shp')


if __name__ == "__main__":
    continental_layers()
