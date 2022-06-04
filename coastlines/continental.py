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

import fiona
import click
import geopandas as gpd
from pathlib import Path

from coastlines.vector import points_on_line
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
    default=[15000, 4000, 1000],
    multiple=True,
    help="The distance (in metres) used to generate coastal "
    "change hotspots summary layers. This controls the spacing "
    "of each summary point, and the radius used to aggregate "
    "rates of change statistics around each point. "
    "The default generates three hotspot layers with radii "
    "15000 m, 4000 m and 1000 m. To specify multiple custom "
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
            f"-single -overwrite_ds -t_srs epsg:6933 "
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
            f"-single -update -t_srs epsg:6933 "
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

        # Clip rates to remove extreme distances, as these are likely due to
        # modelling errors, not true coastal change
        ratesofchange_gdf["rate_time"] = ratesofchange_gdf.rate_time.clip(-250, 250)

        ######################
        # Calculate hotspots #
        ######################

        for i, radius in enumerate(hotspots_radius):

            # Extract hotspot points
            log.info(f"Calculating hotspots at {radius} m")
            hotspots_gdf = points_on_line(
                shorelines_gdf, index=str(baseline_year), distance=radius
            )

            # Create polygon windows
            buffered_gdf = hotspots_gdf[["geometry"]].copy()
            buffered_gdf["geometry"] = buffered_gdf.buffer(radius)

            # Spatial join rate of change points to each polygon, then
            # aggregate/summarise values within each polygon
            hotspot_values = (
                ratesofchange_gdf.sjoin(buffered_gdf, predicate="within")
                .groupby("index_right")["rate_time"]
                .agg([("rate_time", "median"), ("n", "count")])
            )

            # Join aggregated values back to hotspot points
            hotspots_gdf = hotspots_gdf.join(hotspot_values)

            # Add hotspots radius attribute column
            hotspots_gdf["radius_m"] = radius

            # Drop any points with insufficient observations.
            # We can obtain a sensible threshold by dividing the hotspots
            # radius by 30 m along-shore rates of change point distance)
            hotspots_gdf = hotspots_gdf.loc[hotspots_gdf.n > (radius / 30)]

            # Export hotspots to file, incrementing name for each layer
            layer_name = f"hotspots_zoom_{range(0, 10)[i + 1]}"

            try:
                hotspots_gdf.to_file(OUTPUT_FILE, layer=layer_name)
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
