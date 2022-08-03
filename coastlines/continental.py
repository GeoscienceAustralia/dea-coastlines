#!/usr/bin/env python
# coding: utf-8

# This code combines individual datasets into continental DEA
# Coastlines layers:
#
#     * Combines output shorelines and rates of change statistics point
#       vectors into single continental datasets
#     * Aggregates this data to produce moving window hotspot datasets
#       that summarise coastal change at regional and continental scale.


import os
import sys

import fiona
import click
import numpy as np
import geopandas as gpd
from pathlib import Path

from coastlines.vector import points_on_line, change_regress
from coastlines.utils import configure_logging, STYLES_FILE


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
    default=[10000, 5000, 1000],
    multiple=True,
    help="The distance (in metres) used to generate coastal "
    "change hotspots summary layers. This controls the spacing "
    "of each summary point, and the radius used to aggregate "
    "rates of change statistics around each point. "
    "The default generates three hotspot layers with radii "
    "10000 m, 5000 m and 1000 m. To specify multiple custom "
    "radii, repeat this argument, e.g. "
    "`--hotspots_radius 1000 --hotspots_radius 5000`.",
)
@click.option(
    "--baseline_year",
    type=int,
    default=2020,
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

    log = configure_logging("Continental layers and hotspots generation")

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
            f"-single -overwrite_ds -t_srs epsg:3577 "
            f"-nln shorelines_annual"
        )
        log.info("Merging annual shorelines complete")

    else:
        log.info("Not writing shorelines")

    # Combine rates of change stats points into single continental layer
    if ratesofchange:

        os.system(
            f"ogrmerge.py "
            f"-o {OUTPUT_FILE} {ratesofchange_paths} "
            f"-single -update -t_srs epsg:3577 "
            f"-nln rates_of_change"
        )
        log.info("Merging rates of change points complete")

    else:
        log.info("Not writing annual rates of change points")

    #####################
    # Generate hotspots #
    #####################

    # Generate hotspot points that provide regional/continental summary
    # of hotspots of coastal erosion and growth
    if hotspots:

        ###############################
        # Load DEA CoastLines vectors #
        ###############################

        log.info("Generating continental hotspots")

        # Load continental shoreline and rates of change data
        try:

            ratesofchange_gdf = gpd.read_file(OUTPUT_FILE, layer="rates_of_change")
            shorelines_gdf = gpd.read_file(OUTPUT_FILE, layer="shorelines_annual")

        except (fiona.errors.DriverError, ValueError):

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

        # Drop uncertain points from calculation
        ratesofchange_gdf = ratesofchange_gdf.loc[
            ratesofchange_gdf.certainty == "good"
        ].reset_index(drop=True)

        ######################
        # Calculate hotspots #
        ######################

        for i, radius in enumerate(hotspots_radius):

            # Extract hotspot points
            log.info(f"Calculating {radius} m hotspots")
            hotspots_gdf = points_on_line(
                shorelines_gdf, index=str(baseline_year), distance=int(radius / 2)
            )

            # Create polygon windows by buffering points
            buffered_gdf = hotspots_gdf[["geometry"]].copy()
            buffered_gdf["geometry"] = buffered_gdf.buffer(radius)

            # Spatial join rate of change points to each polygon
            hotspot_grouped = (
                ratesofchange_gdf.loc[
                    :, ratesofchange_gdf.columns.str.contains("dist_|geometry")
                ]
                .sjoin(buffered_gdf, predicate="within")
                .groupby("index_right")
            )

            # Aggregate/summarise values by taking median of all points
            # within each buffered polygon
            hotspot_values = hotspot_grouped.median().round(2)

            # Extract year from distance columns (remove "dist_")
            x_years = hotspot_values.columns.str.replace("dist_", "").astype(int)

            # Compute coastal change rates by linearly regressing annual
            # movements vs. time
            rate_out = hotspot_values.apply(
                lambda row: change_regress(
                    y_vals=row.values.astype(float), x_vals=x_years, x_labels=x_years
                ),
                axis=1,
            )

            # Add rates of change back into dataframe
            hotspot_values[
                ["rate_time", "incpt_time", "sig_time", "se_time", "outl_time"]
            ] = rate_out

            # Join aggregated values back to hotspot points after
            # dropping unused columns (regression intercept)
            hotspots_gdf = hotspots_gdf.join(hotspot_values.drop("incpt_time", axis=1))

            # Add hotspots radius attribute column
            hotspots_gdf["radius_m"] = radius

            # Drop any points with insufficient observations.
            # We can obtain a sensible threshold by dividing the hotspots
            # radius by the 30 m along-shore rates of change point distance)
            hotspots_gdf["n"] = hotspot_grouped.size()
            hotspots_gdf = hotspots_gdf.loc[hotspots_gdf.n > (radius / 30)]

            # Export hotspots to file, incrementing name for each layer
            try:

                # Set up schema to optimise file size
                schema_dict = {
                    key: "float:8.2"
                    for key in hotspots_gdf.columns
                    if key != "geometry"
                }
                schema_dict.update(
                    {
                        "sig_time": "float:8.3",
                        "outl_time": "str:80",
                        "radius_m": "int:5",
                        "n": "int:4",
                    }
                )
                col_schema = schema_dict.items()

                # Export file
                layer_name = f"hotspots_zoom_{range(0, 10)[i + 1]}"
                hotspots_gdf.to_file(
                    OUTPUT_FILE,
                    layer=layer_name,
                    schema={"properties": col_schema, "geometry": "Point"},
                )

            except ValueError as e:

                log.exception(f"Failed to generate hotspots with error: {e}")
                sys.exit(1)

        log.info("Writing hotspots complete")

    else:

        log.info("Not writing hotspots...")

    #########################
    # Add GeoPackage styles #
    #########################

    if include_styles:

        log.info("Writing styles in the GeoPackage file")
        styles = gpd.read_file(STYLES_FILE)
        styles.to_file(OUTPUT_FILE, layer="layer_styles")

    else:

        log.info("Not writing styles")


if __name__ == "__main__":
    continental_cli()
