#!/usr/bin/env python
# coding: utf-8

# This code combines individual datasets into continental DEA Coastlines 
# layers:
# 
#     * Combines output coastline and rates of change statistics point 
#       vectors into single continental datasets
#     * Aggregates this data to produce moving window hotspot datasets 
#       that summarise coastal change at regional and continental scale.
#
# Compatability:
#
#     module use /g/data/v10/public/modules/modulefiles
#     module load dea/20200713
#     pip install --user ruptures
#     pip install --user git+https://github.com/mattijn/topojson/


import os
import sys
import geopandas as gpd
import pandas as pd
from rtree import index
from tqdm.auto import tqdm

import deacoastlines_statistics as deacl_stats


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
        poly_ids = [j for j in idx.intersection((point.coords[0]))
                    if point.within(polygons[j])]
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

    matching_points = stats_gdf.iloc[poly_points_dict[key]].copy()

    if len(matching_points.index) > min_n:

        # Set nonsignificant to 0
        matching_points.loc[matching_points.sig_time > 0.01, 'rate_time'] = 0
        matching_points.loc[matching_points.sig_soi > 0.01, 'rate_soi'] = 0

        return pd.Series([matching_points.rate_time.mean(),
                          matching_points.rate_soi.mean(),
                          len(matching_points.index)])

    else:
        return pd.Series([None, None])
    
    
def main(argv=None):
    
    #########
    # Setup #
    #########
    
    if argv is None:

        argv = sys.argv
        print(sys.argv)

    # If no user arguments provided
    if len(argv) < 7:

        str_usage = "You must specify an analysis name"
        print(str_usage)
        sys.exit()
        
    # Set study area and name for analysis
    vector_version = str(argv[1])
    summary_version = str(argv[2])
    threshold = str(argv[3])
    coastlines = bool(argv[4])
    statistics = bool(argv[5])
    summary = argv[6]
    
    #################
    # Merge vectors #
    #################
    
    if coastlines:
        print('Combining annual coastlines')
        os.system(f'ogrmerge.py -o DEACoastlines_annualcoastlines_{summary_version}.shp '
                  f'output_data/*/vectors/shapefiles/contours_*_{vector_version}_'
                  f'mndwi_{threshold}.shp -single -overwrite_ds -t_srs EPSG:3577')
        
    if statistics:
        print('Combining rates of change statistics')
        os.system(f'ogrmerge.py -o DEACoastlines_ratesofchange_{summary_version}.shp '
                  f'output_data/*/vectors/shapefiles/stats_*_{vector_version}_'
                  f'mndwi_{threshold}.shp -single -overwrite_ds -t_srs EPSG:3577')
    
    if summary:

        ###############################
        # Load DEA CoastLines vectors #
        ###############################
        
        print('Generating hotspots')
        stats_gdf = gpd.read_file(f'DEACoastlines_ratesofchange_{summary_version}.shp')
        contours_gdf = gpd.read_file(f'DEACoastlines_annualcoastlines_{summary_version}.shp')

        contours_gdf = (contours_gdf
                        .loc[contours_gdf.geometry.is_valid]
                        .set_index('year'))

        summary_gdf = deacl_stats.points_on_line(contours_gdf, 
                                                 index='2020', 
                                                 distance=summary)

        ####################
        # Generate summary #
        ####################

        # Generate dictionary of polygon IDs and corresponding points
        poly_points_dict = points_in_poly(points=summary_gdf.geometry, 
                                          polygons=stats_gdf.buffer(summary*2))

        # Compute mean and number of obs for each polygon
        summary_gdf[['rate_time', 'rate_soi', 'n']] = summary_gdf.apply(
            lambda row: get_matching_data(row.name, 
                                          stats_gdf,
                                          poly_points_dict,
                                          min_n=summary / 25), axis=1)

        # Export to file
        summary_gdf.to_file(f'DEACoastlines_hotspots_{summary_version}_{summary}.shp')


if __name__ == "__main__":
    main()
