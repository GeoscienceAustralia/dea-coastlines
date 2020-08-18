#!/usr/bin/env python
# coding: utf-8

import os
import sys
import geopandas as gpd
import pandas as pd
from rtree import index
from tqdm.auto import tqdm

import deacoastlines_statistics as deacl_stats


def points_in_poly(points, polygons):

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
    if len(argv) < 6:

        str_usage = "You must specify an analysis name"
        print(str_usage)
        sys.exit()
        
    # Set study area and name for analysis
    output_name = str(argv[1])
    threshold = str(argv[2])
    coastlines = bool(argv[3])
    statistics = bool(argv[4])
    summary = argv[5]
    
    #################
    # Merge vectors #
    #################
    
    if coastlines:
        print('Combining annual coastlines')
        os.system(f'ogrmerge.py -o DEACoastLines_coastlines_{output_name}.shp '
                  f'output_data/*/vectors/shapefiles/contours_*_{output_name}_'
                  f'mndwi_{threshold}.shp -single -overwrite_ds -t_srs EPSG:3577')
        
    if statistics:
        print('Combining rates of change statistics')
        os.system(f'ogrmerge.py -o DEACoastLines_statistics_{output_name}.shp '
                  f'output_data/*/vectors/shapefiles/stats_*_{output_name}_'
                  f'mndwi_{threshold}.shp -single -overwrite_ds -t_srs EPSG:3577')
    
    if summary:

        ###############################
        # Load DEA CoastLines vectors #
        ###############################
        
        print('Generating summary')
        stats_gdf = gpd.read_file(f'DEACoastLines_statistics_{output_name}.shp')
        contours_gdf = gpd.read_file(f'DEACoastLines_coastlines_{output_name}.shp')

        contours_gdf = (contours_gdf
                        .loc[contours_gdf.geometry.is_valid]
                        .set_index('year'))

        summary_gdf = deacl_stats.points_on_line(contours_gdf, 
                                                 index='2019', 
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
        summary_gdf.to_file(f'DEACoastLines_summary_{output_name}_{summary*2}.shp')


if __name__ == "__main__":
    main()