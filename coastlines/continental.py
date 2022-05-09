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

import click
import fiona
import geopandas as gpd
import pandas as pd
from rtree import index
from tqdm.auto import tqdm

from pathlib import Path

from coastlines.vector import points_on_line
from coastlines.utils import configure_logging, STYLES_FILE


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
    for pos, poly in enumerate(tqdm(polygons, desc="Building index")):
        idx.insert(pos, poly.bounds)

    # Iterate through points
    out_dict = {}
    for i, point in enumerate(tqdm(points, desc="Processing points")):
        poly_ids = [
            j for j in idx.intersection((point.coords[0])) if point.within(polygons[j])
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

        return pd.Series([matching_points.rate_time.mean(), len(matching_points.index)])

    else:
        return pd.Series([None, None])


@click.command()
@click.option(
    "--vector_version",
    type=str,
    required=True,
    help="A unique string proving a name that was used "
    "for output vector directories and files. This is used "
    "to identify the tiled annual shoreline and rates of "
    "change layers that will be combined into continental-"
    "scale layers.",
)
@click.option(
    "--continental_version",
    type=str,
    help="A unique string proving a name that will be used "
    "for output continental-scale layers. This allows "
    "multiple versions of continental-scale layers to be "
    "generated from the same input vector data, e.g. for "
    "testing different hotspot of coastal change summary "
    "layers. If not provided, this will default to the "
    'string provided to "--vector_version".',
)
@click.option(
    "--shorelines",
    type=bool,
    default=True,
    help="A boolean indicating whether to combine tiled "
    "annual shorelines layers into a single continental-"
    "scale annual shorelines layer.",
)
@click.option(
    "--ratesofchange",
    type=bool,
    default=True,
    help="A boolean indicating whether to combine tiled "
    "rates of change statistics layers into a single "
    "continental-scale rates of change statistics layer.",
)
@click.option(
    "--hotspots",
    type=bool,
    default=True,
    help="A boolean indicating whether to generate a "
    "continental-scale hotspots of coastal change summary "
    "layer.",
)
@click.option(
    "--hotspots_radius",
    default=10000,
    type=int,
    help="The distance (in metres) used to generate a "
    "coastal change hotspots summary layer. This controls "
    "the spacing of each summary point, and the radius used "
    "to aggregate rates of change statistics around each "
    "point. The default uses a radius of 10000 m.",
)
@click.option(
    "--baseline_year",
    type=str,
    default="2020",
    help="The annual shoreline used to generate the hotspot "
    "summary points. This is typically the most recent "
    "annual shoreline in the dataset.",
)
@click.option(
    "--include-styles/--no-include-styles",
    is_flag=True,
    default=True,
    help="Set this to indicate whether to include styles " "in output Geopackage file.",
)
def continental_cli(
    vector_version,
    continental_version,
    shorelines,
    ratesofchange,
    hotspots,
    hotspots_radius,
    baseline_year,
    include_styles,
):

    #################
    # Merge vectors #
    #################
    log = configure_logging("Coastlines Continental")

    # If no continental version is provided, copy this from vector
    # version
    if continental_version is None:
        continental_version = vector_version
    output_dir = Path(f"data/processed/{continental_version}")
    output_dir.mkdir(exist_ok=True, parents=True)
    log.info(f"Writing data to {output_dir}")

    # Setup input and output file paths
    shoreline_paths = (
        f"data/interim/vector/{vector_version}/*/" f"annualshorelines*.shp"
    )
    ratesofchange_paths = (
        f"data/interim/vector/{vector_version}/*/" f"ratesofchange*.shp"
    )

    OUTPUT_FILE = output_dir / f"coastlines_{continental_version}.gpkg"

    # Combine annual shorelines into a single continental layer
    if shorelines:
        os.system(
            f"ogrmerge.py -o "
            f"{OUTPUT_FILE} {shoreline_paths} "
            f"-single -overwrite_ds -t_srs epsg:6933 "
            f"-nln annual_shorelines"
        )
        log.info("Writing annual shorelines complete")
    else:
        log.info("Not writing shorelines")

    # Combine rates of change stats points into single continental layer
    if ratesofchange:
        os.system(
            f"ogrmerge.py "
            f"-o {OUTPUT_FILE} {ratesofchange_paths} "
            f"-single -update -t_srs epsg:6933 "
            f"-nln rates_of_change"
        )
        log.info("Writing rates of change statistics complete")
    else:
        log.info("Not writing shorelines")

    # Generate hotspot points that provide regional/continental summary
    # of hotspots of coastal erosion and growth
    if hotspots:
        ###############################
        # Load DEA CoastLines vectors #
        ###############################

        log.info("Generating hotspots...")

        # If hotspots is True, set radius from `hotspots_radius`
        hotspots = hotspots_radius

        # Load continental shoreline and rates of change data
        try:
            ratesofchange_gdf = gpd.read_file(OUTPUT_FILE, layer="rates_of_change")
            shorelines_gdf = gpd.read_file(OUTPUT_FILE, layer="annual_shorelines")
        except fiona.errors.DriverError:
            raise FileNotFoundError(
                "Continental-scale annual shoreline and rates of "
                "change layers are required for hotspot generation. "
                "Try re-running this analysis with the following "
                "settings: `--shorelines True --ratesofchange True`."
            )

        # Set year index on coastlines
        shorelines_gdf = shorelines_gdf.loc[shorelines_gdf.geometry.is_valid].set_index(
            "year"
        )

        # Extract hotspot points
        hotspots_gdf = points_on_line(
            shorelines_gdf, index=baseline_year, distance=hotspots
        )

        # Drop low observations (less than 10) from rates
        ratesofchange_gdf = ratesofchange_gdf.loc[ratesofchange_gdf.valid_obs >= 15]
        ratesofchange_gdf = ratesofchange_gdf.reset_index(drop=True)

        # Set nonsignificant rates to 0 m / year
        ratesofchange_gdf.loc[ratesofchange_gdf.sig_time > 0.01, "rate_time"] = 0

        # Clip rates to remove extreme outliers
        ratesofchange_gdf["rate_time"] = ratesofchange_gdf.rate_time.clip(-150, 150)

        #####################
        # Generate hotspots #
        #####################

        # Generate dictionary of polygon IDs and corresponding points
        poly_points_dict = points_in_poly(
            points=hotspots_gdf.geometry,
            polygons=ratesofchange_gdf.buffer(hotspots * 2),
        )

        # Compute mean and number of obs for each polygon
        hotspots_gdf[["rate_time", "n"]] = hotspots_gdf.apply(
            lambda row: get_matching_data(
                row.name, ratesofchange_gdf, poly_points_dict, min_n=hotspots / 30
            ),
            axis=1,
        )

        # Export hotspots to file
        hotspots_gdf.to_file(OUTPUT_FILE, layer="hotspots")
        log.info("Writing hotspots complete")
    else:
        log.info("Not writing hotspots...")

    if include_styles:
        log.info("Writing styles in the geopackage file")
        styles = gpd.read_file(STYLES_FILE)
        styles.to_file(OUTPUT_FILE, layer="layer_styles")
    else:
        log.info("Not writing styles")


if __name__ == "__main__":
    continental_cli()
