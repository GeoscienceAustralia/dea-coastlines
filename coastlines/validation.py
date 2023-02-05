#!/usr/bin/env python
# coding: utf-8

import click
import math
import glob
import re
import os
import os.path
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from pathlib import Path
from io import StringIO
from pyproj import Transformer
from itertools import takewhile
from scipy import stats
import multiprocessing as mp
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from shapely.geometry import box, Point, LineString

import warnings
from shapely.errors import ShapelyDeprecationWarning

warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)


def standardise_source(df):

    # Dictionary containing method values to rename
    remap_dict = {
        "emery": "emery/levelling",
        "levelling": "emery/levelling",
        "dunefill": np.NaN,
        "rtk gps": "gps",
        "photogrammetry": "aerial photogrammetry",
        "stereo photogrammtery": "aerial photogrammetry",
        "ads80": "aerial photogrammetry",
        "photos": "aerial photogrammetry",
        "total station": "total station",
        "total station\t": "total station",
        "laser scanning": "terrestrial laser scanning",
        "satellite": "satellite",
        "gps rtk gps": "gps",
    }

    # Set all values to lower case for easier conversion
    df["source"] = df.source.str.lower()

    # Replace values
    df["source"] = df.source.replace(remap_dict)


def to_vector(
    df, fname="test.shp", x="x", y="y", crs="EPSG:3577", output_crs="EPSG:3577"
):

    # Convert datetimes to strings
    df = df.copy()
    is_datetime = df.dtypes == "datetime64[ns]"
    df.loc[:, is_datetime] = df.loc[:, is_datetime].astype(str)

    # Export to file
    gdf = (
        gpd.GeoDataFrame(
            data=df.loc[:, df.dtypes != "datetime64[ns]"],
            geometry=gpd.points_from_xy(x=df[x], y=df[y]),
            crs=crs,
        )
        .to_crs(output_crs)
        .to_file(fname)
    )

    return gdf


def export_eval(df, output_name, output_crs="EPSG:3577"):

    from shapely.geometry import box, Point, LineString

    # Extract geometries
    val_points = gpd.points_from_xy(x=df.val_x, y=df.val_y)
    deacl_points = gpd.points_from_xy(x=df.deacl_x, y=df.deacl_y)
    df_profiles = df.groupby("id").first()
    profile_lines = df_profiles.apply(
        lambda x: LineString([(x.start_x, x.start_y), (x.end_x, x.end_y)]), axis=1
    )

    # Export validation points
    val_gdf = gpd.GeoDataFrame(data=df, geometry=val_points, crs=output_crs).to_crs(
        "EPSG:4326"
    )
    val_gdf.to_file(f"{output_name}_val.geojson", driver="GeoJSON")

    # Export DEACL points
    deacl_gdf = gpd.GeoDataFrame(data=df, geometry=deacl_points, crs=output_crs).to_crs(
        "EPSG:4326"
    )
    deacl_gdf.to_file(f"{output_name}_deacl.geojson", driver="GeoJSON")

    # Export profiles
    profile_gdf = gpd.GeoDataFrame(
        data=df_profiles, geometry=profile_lines, crs=output_crs
    ).to_crs("EPSG:4326")
    profile_gdf.to_file(f"{output_name}_profiles.geojson", driver="GeoJSON")


def deacl_val_stats(val_dist, deacl_dist, n=None, remove_bias=False):

    np.seterr(all="ignore")

    # Compute difference and bias
    diff_dist = val_dist - deacl_dist
    bias = diff_dist.mean()

    if remove_bias:
        deacl_dist += bias
        diff_dist = val_dist - deacl_dist

    # Compute stats
    if n is None:
        n = len(val_dist)
    else:
        n = sum(n)

    mae = mean_absolute_error(val_dist, deacl_dist)
    rmse = mean_squared_error(val_dist, deacl_dist) ** 0.5

    if n > 1:
        corr = np.corrcoef(x=val_dist, y=deacl_dist)[0][1]
        stdev = diff_dist.std()

    else:
        corr = np.nan
        stdev = np.nan

    return pd.Series(
        {
            "n": n,
            "mae": f"{mae:.2f}",
            "rmse": f"{rmse:.2f}",
            "stdev": f"{stdev:.2f}",
            "corr": f"{corr:.3f}",
            "bias": f"{bias:.2f}",
        }
    ).astype(float)


def rse_tableformat(not_bias_corrected, bias_corrected, groupby="source"):

    # Fix rounding and total observations
    not_bias_corrected["n"] = not_bias_corrected["n"].astype(int)
    not_bias_corrected[["bias", "stdev", "mae", "rmse"]] = not_bias_corrected[
        ["bias", "stdev", "mae", "rmse"]
    ].round(1)
    not_bias_corrected["n"] = not_bias_corrected.groupby(groupby)["n"].sum()

    # Move bias corrected values into brackets
    not_bias_corrected["MAE (m)"] = (
        not_bias_corrected.mae.astype("str")
        + " ("
        + bias_corrected.mae.round(1).astype("str")
        + ")"
    )
    not_bias_corrected["RMSE (m)"] = (
        not_bias_corrected.rmse.astype("str")
        + " ("
        + bias_corrected.rmse.round(1).astype("str")
        + ")"
    )

    # Sort by MAE, rename columns
    not_bias_corrected = (
        not_bias_corrected.sort_values("mae")
        .drop(["mae", "rmse"], axis=1)
        .rename({"stdev": "SD (m)", "corr": "Correlation", "bias": "Bias (m)"}, axis=1)[
            ["n", "Bias (m)", "MAE (m)", "RMSE (m)", "SD (m)", "Correlation"]
        ]
    )

    return not_bias_corrected


def val_slope(profiles_df, intercept_df, datum=0, buffer=25, method="distance"):

    # Join datum dist to full profile dataframe
    profiles_datum_dist = profiles_df.set_index(["id", "date"])[["distance", "z"]].join(
        intercept_df[f"{datum}_dist"]
    )

    if method == "distance":

        # Filter to measurements within distance of datum distance
        beach_data = profiles_datum_dist[
            profiles_datum_dist.distance.between(
                profiles_datum_dist[f"{datum}_dist"] - buffer,
                profiles_datum_dist[f"{datum}_dist"] + buffer,
            )
        ]

    elif method == "height":

        # Filter measurements within height of datum
        beach_data = profiles_datum_dist.loc[
            profiles_datum_dist.z.between(-buffer, buffer)
        ]

    # Calculate slope
    beach_slope = beach_data.groupby(["id", "date"]).apply(
        lambda x: stats.linregress(x=x.distance, y=x.z).slope
    )

    return beach_slope.round(3)


def dms2dd(s):
    # example: s = "0°51'56.29"
    degrees, minutes, seconds = re.split("[°'\"]+", s)
    if float(degrees) > 0:
        dd = float(degrees) + float(minutes) / 60 + float(seconds) / (60 * 60)
    else:
        dd = float(degrees) - float(minutes) / 60 - float(seconds) / (60 * 60)
    return dd


def dist_angle(lon, lat, dist, angle):
    lon_end = lon + dist * np.sin(angle * np.pi / 180)
    lat_end = lat + dist * np.cos(angle * np.pi / 180)
    return pd.Series({"end_y": lat_end, "end_x": lon_end})


def interp_intercept(x, y1, y2, reverse=False):
    """
    Find the intercept of two curves, given by the same x data

    References:
    ----------
    Source: https://stackoverflow.com/a/43551544/2510900
    """

    def intercept(point1, point2, point3, point4):
        """find the intersection between two lines
        the first line is defined by the line between point1 and point2
        the first line is defined by the line between point3 and point4
        each point is an (x,y) tuple.

        So, for example, you can find the intersection between
        intercept((0,0), (1,1), (0,1), (1,0)) = (0.5, 0.5)

        Returns: the intercept, in (x,y) format
        """

        def line(p1, p2):
            A = p1[1] - p2[1]
            B = p2[0] - p1[0]
            C = p1[0] * p2[1] - p2[0] * p1[1]
            return A, B, -C

        def intersection(L1, L2):
            D = L1[0] * L2[1] - L1[1] * L2[0]
            Dx = L1[2] * L2[1] - L1[1] * L2[2]
            Dy = L1[0] * L2[2] - L1[2] * L2[0]

            x = Dx / D
            y = Dy / D
            return x, y

        L1 = line([point1[0], point1[1]], [point2[0], point2[1]])
        L2 = line([point3[0], point3[1]], [point4[0], point4[1]])

        R = intersection(L1, L2)

        return R

    try:

        if isinstance(y2, (int, float)):

            y2 = np.array([y2] * len(x))

        if reverse:

            x = x[::-1]
            y1 = y1[::-1]
            y2 = y2[::-1]

        idx = np.argwhere(np.diff(np.sign(y1 - y2)) != 0)
        xc, yc = intercept(
            (x[idx], y1[idx]),
            ((x[idx + 1], y1[idx + 1])),
            ((x[idx], y2[idx])),
            ((x[idx + 1], y2[idx + 1])),
        )

        return xc[0][0]

    except:

        return np.nan


def dist_along_transect(dist, start_x, start_y, end_x, end_y):

    transect_line = LineString([(start_x, start_y), (end_x, end_y)])
    distance_coords = transect_line.interpolate(dist).coords.xy
    return [coord[0] for coord in distance_coords]


def waterline_intercept(
    x, dist_col="distance", x_col="x", y_col="y", z_col="z", z_val=0, debug=False
):

    # Extract distance and coordinates of where the z_val first
    # intersects with the profile line
    dist_int = interp_intercept(x[dist_col].values, x[z_col].values, z_val)
    x_int = interp_intercept(x[x_col].values, x[z_col].values, z_val)
    y_int = interp_intercept(x[y_col].values, x[z_col].values, z_val)

    # Identify last distance where the z_value intersects the profile
    rev_int = interp_intercept(x[dist_col].values, x[z_col].values, z_val, reverse=True)

    # If first and last intersects are the identical, return data.
    # If not, the comparison is invalid (i.e. NaN)
    if dist_int == rev_int:
        if debug:
            print("Single intersection found")
        return pd.Series(
            {f"{z_val}_dist": dist_int, f"{z_val}_x": x_int, f"{z_val}_y": y_int}
        )
    else:
        if debug:
            print("Multiple intersections returned")
        return pd.Series(
            {f"{z_val}_dist": np.NaN, f"{z_val}_x": np.NaN, f"{z_val}_y": np.NaN}
        )


def reproj_crs(in_data, in_crs, x="x", y="y", out_crs="EPSG:3577"):

    # Reset index to allow merging new data with original data
    in_data = in_data.reset_index(drop=True)

    # Reproject coords to Albers and create geodataframe
    trans = Transformer.from_crs(in_crs, out_crs, always_xy=True)
    coords = trans.transform(in_data[x].values, in_data[y].values)
    in_data[["x", "y"]] = pd.DataFrame(zip(*coords))

    return in_data


def profiles_from_dist(
    profiles_df, id_col="id", dist_col="distance", x_col="x", y_col="y"
):

    # Compute origin points for each profile
    min_ids = profiles_df.groupby(id_col)[dist_col].idxmin()
    start_xy = profiles_df.loc[min_ids, [id_col, x_col, y_col]]
    start_xy = start_xy.rename(
        {x_col: f"start_{x_col}", y_col: f"start_{y_col}"}, axis=1
    )

    # Compute end points for each profile
    max_ids = profiles_df.groupby(id_col)[dist_col].idxmax()
    end_xy = profiles_df.loc[max_ids, [x_col, y_col]]

    # Add end coords into same dataframe
    start_xy = start_xy.reset_index(drop=True)
    end_xy = end_xy.reset_index(drop=True)
    start_xy[[f"end_{x_col}", f"end_{y_col}"]] = end_xy

    return start_xy


def perpendicular_line(input_line, length):

    # Generate parallel lines either side of input line
    left = input_line.parallel_offset(length / 2.0, "left")
    right = input_line.parallel_offset(length / 2.0, "right")

    # Create new line between centroids of parallel line.
    # This should be perpendicular to the original line
    return LineString([left.centroid, right.centroid])


def generate_transects(line_geom, length=400, interval=200, buffer=20):

    # Create tangent line at equal intervals along line geom
    interval_dists = np.arange(buffer, line_geom.length, interval)
    tangent_geom = [
        LineString(
            [line_geom.interpolate(dist - buffer), line_geom.interpolate(dist + buffer)]
        )
        for dist in interval_dists
    ]

    # Convert to geoseries and remove erroneous lines by length
    tangent_gs = gpd.GeoSeries(tangent_geom)
    tangent_gs = tangent_gs.loc[tangent_gs.length.round(1) <= buffer * 2]

    # Compute perpendicular lines
    return tangent_gs.apply(perpendicular_line, length=length)


def coastal_transects(
    bbox,
    name,
    interval=200,
    transect_length=400,
    simplify_length=200,
    transect_buffer=20,
    output_crs="EPSG:3577",
    coastline="../input_data/Smartline.gdb",
    land_poly="/g/data/r78/rt1527/shapefiles/australia/australia/cstauscd_r.shp",
):

    # Load smartline
    coastline_gdf = gpd.read_file(coastline, bbox=bbox).to_crs(output_crs)
    coastline_geom = coastline_gdf.geometry.unary_union.simplify(simplify_length)

    # Load Australian land water polygon
    land_gdf = gpd.read_file(land_poly, bbox=bbox).to_crs(output_crs)
    land_gdf = land_gdf.loc[land_gdf.FEAT_CODE.isin(["mainland", "island"])]
    land_geom = gpd.overlay(df1=land_gdf, df2=bbox).unary_union

    # Extract transects along line
    geoms = generate_transects(
        coastline_geom,
        length=transect_length,
        interval=interval,
        buffer=transect_buffer,
    )

    # Test if end points of transects fall in water or land
    p1 = gpd.GeoSeries([Point(i.coords[0]) for i in geoms])
    p2 = gpd.GeoSeries([Point(i.coords[1]) for i in geoms])
    p1_within_land = p1.within(land_geom)
    p2_within_land = p2.within(land_geom)

    # Create geodataframe, remove invalid land-land/water-water transects
    transect_gdf = gpd.GeoDataFrame(
        data={"p1": p1_within_land, "p2": p2_within_land},
        geometry=geoms.values,
        crs=output_crs,
    )
    transect_gdf = transect_gdf[~(transect_gdf.p1 == transect_gdf.p2)]

    # Reverse transects so all point away from land
    transect_gdf["geometry"] = transect_gdf.apply(
        lambda i: LineString([i.geometry.coords[1], i.geometry.coords[0]])
        if i.p1 < i.p2
        else i.geometry,
        axis=1,
    )

    # Export to file
    transect_gdf[["geometry"]].to_file(
        f"input_data/coastal_transects_{name}.geojson", driver="GeoJSON"
    )


def coastal_transects_parallel(
    regions_gdf,
    interval=200,
    transect_length=400,
    simplify_length=200,
    transect_buffer=20,
    overwrite=False,
    output_path="input_data/combined_transects_wadot.geojson",
):

    if not os.path.exists(output_path) or overwrite:

        if os.path.exists(output_path):
            print("Removing existing file")
            os.remove(output_path)

        # Generate transects for each region
        print("Generating transects")
        with mp.Pool(mp.cpu_count()) as pool:
            for i, _ in regions_gdf.iterrows():
                name = str(i).replace(" ", "").replace("/", "").lower()
                pool.apply_async(
                    coastal_transects,
                    [
                        regions_gdf.loc[[i]],
                        name,
                        interval,
                        transect_length,
                        simplify_length,
                        transect_buffer,
                    ],
                )
            pool.close()
            pool.join()

        # Load regional transects and combine into a single file
        print("Combining data")
        transect_list = glob.glob("input_data/coastal_transects_*.geojson")
        gdf = pd.concat(
            [gpd.read_file(shp, ignore_index=True) for shp in transect_list]
        )
        gdf = gdf.reset_index(drop=True)
        gdf["profile"] = gdf.index.astype(str)
        gdf.to_file(output_path, driver="GeoJSON")

        # Clean files
        [os.remove(f) for f in transect_list]


def preprocess_wadot(
    compartment,
    overwrite=True,
    fname="input_data/wadot/Coastline_Movements_20190819.gdb",
):

    beach = str(compartment.index.item())
    fname_out = f"output_data/wadot_{beach}.csv"
    print(f"Processing {beach:<80}", end="\r")

    # Test if file exists
    if not os.path.exists(fname_out) or overwrite:

        # Read file and filter to AHD 0 shorelines
        val_gdf = gpd.read_file(fname, bbox=compartment).to_crs("EPSG:3577")
        val_gdf = gpd.clip(gdf=val_gdf, mask=compartment, keep_geom_type=True)
        val_gdf = val_gdf[(val_gdf.TYPE == "AHD 0m") | (val_gdf.TYPE == "AHD 0m ")]

        # Filter to post 1987 shorelines and set index to year
        val_gdf = val_gdf[val_gdf.PHOTO_YEAR > 1987]
        val_gdf = val_gdf.set_index("PHOTO_YEAR")

        # If no data is returned, skip this iteration
        if len(val_gdf.index) == 0:
            print(f"Failed: {beach:<80}", end="\r")
            return None

        ######################
        # Generate transects #
        ######################

        transect_gdf = gpd.read_file(
            "input_data/combined_transects_wadot.geojson", bbox=compartment
        )
        transect_gdf = gpd.clip(gdf=transect_gdf, mask=compartment, keep_geom_type=True)

        ################################
        # Identify 0 MSL intersections #
        ################################

        output_list = []

        # Select one year
        for year in val_gdf.index.unique().sort_values():

            # Extract validation contour
            print(f"Processing {beach} {year:<80}", end="\r")
            val_contour = val_gdf.loc[[year]].geometry.unary_union

            # Copy transect data, and find intersects
            # between transects and contour
            intersect_gdf = transect_gdf.copy()
            intersect_gdf["val_point"] = transect_gdf.intersection(val_contour)
            to_keep = gpd.GeoSeries(intersect_gdf["val_point"]).geom_type == "Point"
            intersect_gdf = intersect_gdf.loc[to_keep]

            # If no data is returned, skip this iteration
            if len(intersect_gdf.index) == 0:
                print(f"Failed: {beach} {year:<80}", end="\r")
                continue

            # Add generic metadata
            intersect_gdf["date"] = pd.to_datetime(str(year))
            intersect_gdf["beach"] = beach
            intersect_gdf["section"] = "all"
            intersect_gdf["source"] = "aerial photogrammetry"
            intersect_gdf["name"] = "wadot"
            intersect_gdf["id"] = (
                intersect_gdf.beach
                + "_"
                + intersect_gdf.section
                + "_"
                + intersect_gdf.profile
            )

            # Add measurement metadata
            intersect_gdf[["start_x", "start_y"]] = intersect_gdf.apply(
                lambda x: pd.Series(x.geometry.coords[0]), axis=1
            )
            intersect_gdf[["end_x", "end_y"]] = intersect_gdf.apply(
                lambda x: pd.Series(x.geometry.coords[1]), axis=1
            )
            intersect_gdf["0_dist"] = intersect_gdf.apply(
                lambda x: Point(x.start_x, x.start_y).distance(x["val_point"]), axis=1
            )
            intersect_gdf[["0_x", "0_y"]] = intersect_gdf.apply(
                lambda x: pd.Series(x.val_point.coords[0]), axis=1
            )

            # Add empty slope var (not possible to compute without profile data)
            intersect_gdf["slope"] = np.nan

            # Keep required columns
            intersect_gdf = intersect_gdf[
                [
                    "id",
                    "date",
                    "beach",
                    "section",
                    "profile",
                    "name",
                    "source",
                    "slope",
                    "start_x",
                    "start_y",
                    "end_x",
                    "end_y",
                    "0_dist",
                    "0_x",
                    "0_y",
                ]
            ]

            # Append to file
            output_list.append(intersect_gdf)

        # Combine all year data and export to file
        if len(output_list) > 0:
            shoreline_df = pd.concat(output_list)
            shoreline_df.to_csv(fname_out, index=False)

    else:
        print(f"Skipping {beach:<80}", end="\r")


def preprocess_dasilva2021(
    fname="input_data/dasilva2021/dasilva_etal_2021_shorelines.shp",
):

    beach = "dasilva2021"
    print(f"Processing {beach:<80}", end="\r")

    # Read file and filter to AHD 0 shorelines
    fname = "input_data/dasilva2021/dasilva_etal_2021_shorelines.shp"
    val_gdf = gpd.read_file(fname).to_crs("EPSG:3577")
    val_gdf = val_gdf.loc[val_gdf.Year_ > 1987]
    val_gdf["Year_"] = val_gdf.Year_.astype(str)
    val_gdf = val_gdf.set_index("Year_")

    # If no data is returned, skip this iteration
    if len(val_gdf.index) == 0:
        print(f"Failed: {beach:<80}", end="\r")
        return None

    ######################
    # Generate transects #
    ######################

    transect_gdf = gpd.read_file(
        "input_data/dasilva2021/dasilva_etal_2021_retransects.shp"
    ).to_crs("EPSG:3577")[["TransectID", "Direction", "order", "geometry"]]
    transect_gdf.columns = ["profile", "section", "order", "geometry"]
    transect_gdf = transect_gdf.sort_values("order").set_index("order")
    transect_gdf["profile"] = transect_gdf.profile.astype(str)

    ################################
    # Identify 0 MSL intersections #
    ################################

    output_list = []

    # Select one year
    for year in val_gdf.index.unique().sort_values():

        # Extract validation contour
        print(f"Processing {beach} {year:<80}", end="\r")
        val_contour = val_gdf.loc[[year]].geometry.unary_union

        # Copy transect data, and find intersects
        # between transects and contour
        intersect_gdf = transect_gdf.copy()
        intersect_gdf["val_point"] = transect_gdf.intersection(val_contour)
        to_keep = gpd.GeoSeries(intersect_gdf["val_point"]).geom_type == "Point"
        intersect_gdf = intersect_gdf.loc[to_keep]

        # If no data is returned, skip this iteration
        if len(intersect_gdf.index) == 0:
            print(f"Failed: {beach} {year:<80}", end="\r")
            continue

        # Add generic metadata
        intersect_gdf["date"] = pd.to_datetime(str(year))
        intersect_gdf["beach"] = beach
        intersect_gdf["source"] = "satellite"
        intersect_gdf["name"] = "dasilva2021"
        intersect_gdf["id"] = (
            intersect_gdf.beach
            + "_"
            + intersect_gdf.section
            + "_"
            + intersect_gdf.profile
        )

        # Add measurement metadata
        intersect_gdf[["start_x", "start_y"]] = intersect_gdf.apply(
            lambda x: pd.Series(x.geometry.coords[0]), axis=1
        )
        intersect_gdf[["end_x", "end_y"]] = intersect_gdf.apply(
            lambda x: pd.Series(x.geometry.coords[1]), axis=1
        )
        intersect_gdf["0_dist"] = intersect_gdf.apply(
            lambda x: Point(x.start_x, x.start_y).distance(x["val_point"]), axis=1
        )
        intersect_gdf[["0_x", "0_y"]] = intersect_gdf.apply(
            lambda x: pd.Series(x.val_point.coords[0][0:2]), axis=1
        )

        # Add empty slope var (not possible to compute without profile data)
        intersect_gdf["slope"] = np.nan

        # Keep required columns
        intersect_gdf = intersect_gdf[
            [
                "id",
                "date",
                "beach",
                "section",
                "profile",
                "name",
                "source",
                "slope",
                "start_x",
                "start_y",
                "end_x",
                "end_y",
                "0_dist",
                "0_x",
                "0_y",
            ]
        ]

        # Append to file
        output_list.append(intersect_gdf)

    # Combine all year data and export to file
    if len(output_list) > 0:
        shoreline_df = pd.concat(output_list)
        shoreline_df.to_csv(f"output_data/{beach}.csv", index=False)


def preprocess_stirling(fname_out, datum=0):

    # List containing files to import and all params to extract data
    survey_xl = [
        {
            "fname": "input_data/stirling/2015 05 28 - From Stirling - Coastal Profiles 2014-2015 April-Feb with updated reef#2.xlsm",
            "skiprows": 2,
            "skipcols": 5,
            "nrows": 100,
            "meta_skiprows": 0,
            "meta_nrows": 1,
            "meta_usecols": [6, 7],
        },
        {
            "fname": "input_data/stirling/Coastal Profiles 2013-2014 JUL-MAY#2.xlsx",
            "skiprows": 2,
            "skipcols": 5,
            "nrows": 100,
            "meta_skiprows": 0,
            "meta_nrows": 1,
            "meta_usecols": [6, 7],
        },
        {
            "fname": "input_data/stirling/COASTAL PROFILES 2013 JAN - JUNE#2.xls",
            "skiprows": 3,
            "skipcols": 0,
            "nrows": 40,
            "meta_skiprows": 1,
            "meta_nrows": 2,
            "meta_usecols": [1, 2],
        },
        {
            "fname": "input_data/stirling/COASTAL PROFILES 2012 JUN - DEC#2.xls",
            "skiprows": 3,
            "skipcols": 0,
            "nrows": 40,
            "meta_skiprows": 1,
            "meta_nrows": 2,
            "meta_usecols": [1, 2],
        },
        {
            "fname": "input_data/stirling/COASTAL PROFILES 2011-2012 NOV - MAY#2.xls",
            "skiprows": 3,
            "skipcols": 0,
            "nrows": 40,
            "meta_skiprows": 1,
            "meta_nrows": 2,
            "meta_usecols": [1, 2],
        },
    ]

    # List to contain processed profile data
    output = []

    # For each survey excel file in the list above:
    for survey in survey_xl:

        # Load profile start metadata
        all_meta = pd.read_excel(
            survey["fname"],
            sheet_name=None,
            nrows=survey["meta_nrows"],
            skiprows=survey["meta_skiprows"],
            usecols=survey["meta_usecols"],
            header=None,
            on_demand=True,
        )

        # Load data
        all_sheets = pd.read_excel(
            survey["fname"],
            sheet_name=None,
            skiprows=survey["skiprows"],
            nrows=survey["nrows"],
            parse_dates=False,
            usecols=lambda x: "Unnamed" not in str(x),
        )

        # Iterate through each profile in survey data
        for profile_id in np.arange(1, 20).astype("str"):

            # Extract profile start metadata and profile data
            start_x, start_y = all_meta[profile_id].values[0]
            sheet = all_sheets[profile_id].iloc[:, survey["skipcols"] :]

            # First set all column names to lower case strings
            sheet.columns = sheet.columns.astype(str).str.slice(0, 10).str.lower()

            # Drop note columns and distance/angle offset
            sheet = sheet.loc[:, ~sheet.columns.str.contains("note|notes")]
            sheet = sheet.drop(["dist", "angle dd"], axis=1, errors="ignore")

            # Expand date column values into rows for each sampling event
            sheet.loc[:, sheet.columns[::4]] = sheet.columns[::4]

            # Number date columns incrementally to match other fields
            start_num = 1 if survey["skipcols"] > 0 else 0
            rename_dict = {
                name: f"date.{i + start_num}"
                for i, name in enumerate(sheet.columns[::4])
            }
            sheet = sheet.rename(rename_dict, axis=1).reset_index()
            sheet = sheet.rename({"x": "x.0", "y": "y.0", "z": "z.0"}, axis=1)

            # Reshape data into long format
            profile_df = pd.wide_to_long(
                sheet, stubnames=["date", "x", "y", "z"], i="index", j="dropme", sep="."
            ).reset_index(drop=True)

            # Set datetimes
            profile_df["date"] = pd.to_datetime(
                profile_df.date, errors="coerce", dayfirst=True
            )

            # Add profile metadata
            profile_df["beach"] = "stirling"
            profile_df["section"] = "all"
            profile_df["profile"] = profile_id
            profile_df["name"] = "stirling"
            profile_df["source"] = "gps"
            profile_df["start_x"] = start_x
            profile_df["start_y"] = start_y
            profile_df["id"] = (
                profile_df.beach + "_" + profile_df.section + "_" + profile_df.profile
            )

            # Add results to list
            output.append(profile_df.dropna())

    # Combine all survey and profile data
    profiles_df = pd.concat(output)

    # Reproject Perth Coastal Grid coordinates into Australian Albers
    pcg_crs = (
        "+proj=tmerc +lat_0=0 +lon_0=115.8166666666667 "
        "+k=0.9999990600000001 +x_0=50000 +y_0=3800000 "
        "+ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
    )
    trans = Transformer.from_crs(pcg_crs, "EPSG:3577", always_xy=True)
    profiles_df["x"], profiles_df["y"] = trans.transform(
        profiles_df.y.values, profiles_df.x.values
    )
    profiles_df["start_x"], profiles_df["start_y"] = trans.transform(
        profiles_df.start_y.values, profiles_df.start_x.values
    )

    # Calculate per-point distance from start of profile
    profiles_df["distance"] = profiles_df.apply(
        lambda x: Point(x.start_x, x.start_y).distance(Point(x.x, x.y)), axis=1
    )

    # Identify end of profiles by max distance from start, and merge back
    max_dist = (
        profiles_df.sort_values("distance", ascending=False)
        .groupby("id")["x", "y"]
        .first()
        .rename({"x": "end_x", "y": "end_y"}, axis=1)
    )
    profiles_df = profiles_df.merge(max_dist, on="id")

    # Find location and distance to water for datum height (e.g. 0 m AHD)
    intercept_df = (
        profiles_df.groupby(["id", "date"])
        .apply(waterline_intercept, z_val=datum)
        .dropna()
    )

    # Join into dataframe
    shoreline_dist = intercept_df.join(profiles_df.groupby(["id", "date"]).first())

    # Keep required columns
    shoreline_dist = shoreline_dist[
        [
            "beach",
            "section",
            "profile",
            "name",
            "source",
            "start_x",
            "start_y",
            "end_x",
            "end_y",
            f"{datum}_dist",
            f"{datum}_x",
            f"{datum}_y",
        ]
    ]

    # Export to file
    shoreline_dist.to_csv(fname_out)


def preprocess_vicdeakin(fname, datum=0):

    # Dictionary to map correct CRSs to locations
    crs_dict = {
        "apo": "epsg:32754",
        "cow": "epsg:32755",
        "inv": "epsg:32755",
        "leo": "epsg:32755",
        "mar": "epsg:32754",
        "pfa": "epsg:32754",
        "por": "epsg:32755",
        "prd": "epsg:32755",
        "sea": "epsg:32755",
        "wbl": "epsg:32754",
    }

    # Read data
    profiles_df = pd.read_csv(fname, parse_dates=["survey_date"]).dropna()

    # Restrict to pre-2019
    profiles_df = profiles_df.loc[profiles_df.survey_date.dt.year < 2020]
    profiles_df = profiles_df.reset_index(drop=True)

    # Remove invalid profiles
    invalid = (profiles_df.location == "leo") & (profiles_df.tr_id == 94)
    profiles_df = profiles_df.loc[~invalid].reset_index(drop=True)

    # Extract coordinates
    coords = profiles_df.coordinates.str.findall(r"\d+\.\d+")
    profiles_df[["x", "y"]] = pd.DataFrame(coords.values.tolist(), dtype=np.float32)

    # Add CRS and convert to Albers
    profiles_df["crs"] = profiles_df.location.apply(lambda x: crs_dict[x])
    profiles_df = (
        profiles_df.groupby("crs", as_index=False)
        .apply(lambda x: reproj_crs(x, in_crs=x.crs.iloc[0]))
        .drop("crs", axis=1)
    )
    profiles_df = profiles_df.reset_index(drop=True)

    # Convert columns to strings and add unique ID column
    profiles_df = profiles_df.rename(
        {
            "location": "beach",
            "tr_id": "profile",
            "survey_date": "date",
            "z": "z_dirty",
            "z_clean": "z",
        },
        axis=1,
    )
    profiles_df["profile"] = profiles_df["profile"].astype(str)
    profiles_df["section"] = "all"
    profiles_df["source"] = "drone photogrammetry"
    profiles_df["name"] = "vicdeakin"
    profiles_df["id"] = (
        profiles_df.beach + "_" + profiles_df.section + "_" + profiles_df.profile
    )

    # Reverse profile distances by subtracting max distance from each
    prof_max = profiles_df.groupby("id")["distance"].transform("max")
    profiles_df["distance"] = (profiles_df["distance"] - prof_max).abs()

    # Compute origin and end points for each profile and merge into data
    start_end_xy = profiles_from_dist(profiles_df)
    profiles_df = pd.merge(left=profiles_df, right=start_end_xy)

    # Export each beach
    for beach_name, beach in profiles_df.groupby("beach"):

        # Create output file name
        fname_out = f"output_data/vicdeakin_{beach_name}.csv"
        print(f"Processing {fname_out:<80}", end="\r")

        # Find location and distance to water for datum height (0 m AHD)
        intercept_df = (
            beach.groupby(["id", "date"])
            .apply(waterline_intercept, z_val=datum)
            .dropna()
        )

        # If the output contains data
        if len(intercept_df.index) > 0:

            # Join into dataframe
            shoreline_dist = intercept_df.join(beach.groupby(["id", "date"]).first())

            # Compute validation slope and join into dataframe
            slope = val_slope(beach, intercept_df, datum=datum)
            shoreline_dist = shoreline_dist.join(slope.rename("slope"))

            # Keep required columns
            shoreline_dist = shoreline_dist[
                [
                    "beach",
                    "section",
                    "profile",
                    "name",
                    "source",
                    "slope",
                    "start_x",
                    "start_y",
                    "end_x",
                    "end_y",
                    f"{datum}_dist",
                    f"{datum}_x",
                    f"{datum}_y",
                ]
            ]

        # Export to file
        shoreline_dist.to_csv(fname_out)


def preprocess_nswbpd(fname, datum=0, overwrite=False):

    # Get output filename
    name = Path(fname).stem.split("_")[-1].lower().replace(" ", "")
    fname_out = f"output_data/nswbpd_{name}.csv"

    # Test if file exists
    if not os.path.exists(fname_out) or overwrite:

        # Read in data
        print(f"Processing {fname_out:<80}", end="\r")
        profiles_df = pd.read_csv(
            fname, skiprows=5, dtype={"Block": str, "Profile": str}
        )
        profiles_df["Year/Date"] = pd.to_datetime(
            profiles_df["Year/Date"], dayfirst=True, errors="coerce"
        )

        # Convert columns to strings and add unique ID column
        profiles_df["Beach"] = profiles_df["Beach"].str.lower().str.replace(" ", "")
        profiles_df["Block"] = profiles_df["Block"].str.lower()
        profiles_df["Profile"] = profiles_df["Profile"].astype(str).str.lower()
        profiles_df["id"] = (
            profiles_df.Beach + "_" + profiles_df.Block + "_" + profiles_df.Profile
        )
        profiles_df["name"] = "nswbpd"

        # Rename columns
        profiles_df.columns = [
            "beach",
            "section",
            "profile",
            "date",
            "distance",
            "z",
            "x",
            "y",
            "source",
            "id",
            "name",
        ]

        # Reproject coords to Albers
        trans = Transformer.from_crs("EPSG:32756", "EPSG:3577", always_xy=True)
        profiles_df["x"], profiles_df["y"] = trans.transform(
            profiles_df.x.values, profiles_df.y.values
        )

        # Restrict to post 1987
        profiles_df = profiles_df[profiles_df["date"] > "1987"]

        # Compute origin and end points for each profile and merge into data
        start_end_xy = profiles_from_dist(profiles_df)
        profiles_df = pd.merge(left=profiles_df, right=start_end_xy)

        # Drop profiles that have been assigned incorrect profile IDs.
        # To do this, we use a correlation test to determine whether x
        # and y coordinates within each individual profiles fall along a
        # straight line. If a profile has a low correlation (e.g. less
        # than 99.9), it is likely that multiple profile lines have been
        # incorrectly labelled with a single profile ID.
        valid_profiles = lambda x: x[["x", "y"]].corr().abs().iloc[0, 1] > 0.99
        drop = (~profiles_df.groupby("id").apply(valid_profiles)).sum()
        profiles_df = profiles_df.groupby("id").filter(valid_profiles)
        if drop.sum() > 0:
            print(f"\nDropping invalid profiles: {drop:<80}")

        # If profile data remains
        if len(profiles_df.index) > 0:

            # Restrict profiles to data that falls ocean-ward of the top of
            # the foredune (the highest point in the profile) to remove
            # spurious validation points, e.g. due to a non-shoreline lagoon
            # at the back of the profile
            foredune_dist = (
                profiles_df.groupby(["id", "date"])
                .apply(lambda x: x.distance.loc[x.z.idxmax()])
                .reset_index(name="foredune_dist")
            )
            profiles_df = pd.merge(left=profiles_df, right=foredune_dist)
            profiles_df = profiles_df.loc[
                (profiles_df.distance >= profiles_df.foredune_dist)
            ]

            # Find location and distance to water for datum height (e.g. 0 m AHD)
            intercept_df = (
                profiles_df.groupby(["id", "date"])
                .apply(waterline_intercept, z_val=datum)
                .dropna()
            )

            # If any datum intercepts are found
            if len(intercept_df.index) > 0:

                # Join into dataframe
                shoreline_dist = intercept_df.join(
                    profiles_df.groupby(["id", "date"]).agg(
                        lambda x: pd.Series.mode(x).iloc[0]
                    )
                )

                # Compute validation slope and join into dataframe
                slope = val_slope(profiles_df, intercept_df, datum=datum)
                shoreline_dist = shoreline_dist.join(slope.rename("slope"))

                # Keep required columns
                shoreline_dist = shoreline_dist[
                    [
                        "beach",
                        "section",
                        "profile",
                        "name",
                        "source",
                        "foredune_dist",
                        "slope",
                        "start_x",
                        "start_y",
                        "end_x",
                        "end_y",
                        f"{datum}_dist",
                        f"{datum}_x",
                        f"{datum}_y",
                    ]
                ]

                # Standardise source column
                standardise_source(shoreline_dist)

                # Export to file
                shoreline_dist.to_csv(fname_out)

    else:
        print(f"Skipping {fname:<80}", end="\r")


def preprocess_narrabeen(
    fname, fname_out="output_data/wrl_narrabeen.csv", datum=0, overwrite=False
):

    # Test if file exists
    if not os.path.exists(fname_out) or overwrite:

        #################
        # Location data #
        #################

        # Import data and parse DMS to DD
        print(f"Processing {fname_out:<80}", end="\r")
        data = (
            "PF1 -33°42'20.65 151°18'16.30 118.42\n"
            "PF2 -33°42'33.45 151°18'10.33 113.36\n"
            "PF4 -33°43'01.55 151°17'58.84 100.26\n"
            "PF6 -33°43'29.81 151°17'58.65 83.65\n"
            "PF8 -33°43'55.94 151°18'06.47 60.48"
        )
        coords = pd.read_csv(
            StringIO(data), sep=" ", names=["profile", "y", "x", "angle"]
        )
        coords["x"] = [dms2dd(i) for i in coords.x]
        coords["y"] = [dms2dd(i) for i in coords.y]

        # Extend survey lines out from start coordinates using supplied angle
        coords_end = coords.apply(
            lambda x: dist_angle(x.x, x.y, 0.002, x.angle), axis=1
        )
        coords = pd.concat([coords, coords_end], axis=1).drop("angle", axis=1)

        # Rename initial x and y values to indicate profile starting coords
        coords = coords.rename({"y": "start_y", "x": "start_x"}, axis=1)

        # Reproject coords to Albers and create geodataframe
        trans = Transformer.from_crs("EPSG:4326", "EPSG:3577", always_xy=True)
        coords["start_x"], coords["start_y"] = trans.transform(
            coords.start_x.values, coords.start_y.values
        )
        coords["end_x"], coords["end_y"] = trans.transform(
            coords.end_x.values, coords.end_y.values
        )

        # Add ID column
        coords["profile"] = coords["profile"].astype(str).str.lower()
        coords["beach"] = "narrabeen"
        coords["section"] = "all"
        coords["name"] = "wrl"
        coords["id"] = coords.beach + "_" + coords.section + "_" + coords.profile

        ###############
        # Survey data #
        ###############

        # Import data
        profiles_df = pd.read_csv(
            fname,
            usecols=[1, 2, 3, 4, 5],
            skiprows=1,
            parse_dates=["date"],
            names=["profile", "date", "distance", "z", "source"],
        )

        # Restrict to post 1987
        profiles_df = profiles_df[(profiles_df.date.dt.year > 1987)]

        # Merge profile coordinate data into transect data
        profiles_df["profile"] = profiles_df["profile"].astype(str).str.lower()
        profiles_df = profiles_df.merge(coords, on="profile")

        # Add coordinates at every supplied distance along transects
        profiles_df[["x", "y"]] = profiles_df.apply(
            lambda x: pd.Series(
                dist_along_transect(x.distance, x.start_x, x.start_y, x.end_x, x.end_y)
            ),
            axis=1,
        )

        # Find location and distance to water for datum height (e.g. 0 m AHD)
        intercept_df = (
            profiles_df.groupby(["id", "date"])
            .apply(waterline_intercept, z_val=datum)
            .dropna()
        )

        # If the output contains data
        if len(intercept_df.index) > 0:

            # Join into dataframe
            shoreline_dist = intercept_df.join(
                profiles_df.groupby(["id", "date"]).agg(
                    lambda x: pd.Series.mode(x).iloc[0]
                )
            )

            # Compute validation slope and join into dataframe
            slope = val_slope(profiles_df, intercept_df, datum=datum)
            shoreline_dist = shoreline_dist.join(slope.rename("slope"))

            # Keep required columns
            shoreline_dist = shoreline_dist[
                [
                    "beach",
                    "section",
                    "profile",
                    "name",
                    "source",
                    "slope",
                    "start_x",
                    "start_y",
                    "end_x",
                    "end_y",
                    f"{datum}_dist",
                    f"{datum}_x",
                    f"{datum}_y",
                ]
            ]

            # Standardise source column
            standardise_source(shoreline_dist)

            # Export to file
            shoreline_dist.to_csv(fname_out)

    else:
        print(f"Skipping {fname:<80}", end="\r")


def preprocess_cgc(site, datum=0, overwrite=True):

    # Standardise beach name from site name
    beach = site.replace("NO*TH KIRRA", "NORTH KIRRA").lower()
    beach = beach.replace(" ", "").lower()
    fname_out = f"output_data/cgc_{beach}.csv"
    print(f"Processing {fname_out:<80}", end="\r")

    # Test if file exists
    if not os.path.exists(fname_out) or overwrite:

        # List of profile datasets to iterate through
        profile_list = glob.glob(f"input_data/cityofgoldcoast/{site}*.txt")

        # Output list to hold data
        site_profiles = []

        # For each profile, import and standardise data then add to list
        for profile_i in profile_list:

            # Identify unique field values from file string
            profile_string = os.path.basename(profile_i)
            date = profile_string.split(" - (")[1][-14:-4]
            section_profile = profile_string.split(" - (")[0].split(" - ")[1]
            section = section_profile.split(" ")[0]
            profile = "".join(section_profile.split(" ")[1:])

            # Fix missing section or profile info
            if section and not profile:
                section, profile = "na", section
            elif not section and not profile:
                section, profile = "na", "na"

            # Set location metadata and ID
            profile_df = pd.read_csv(
                profile_i,
                usecols=[1, 2, 3],
                delim_whitespace=True,
                names=["x", "y", "z"],
            )
            profile_df["date"] = pd.to_datetime(date)
            profile_df["source"] = "hydrographic survey"
            profile_df["name"] = "cgc"
            profile_df["profile"] = profile.lower()
            profile_df["section"] = section.lower()
            profile_df["beach"] = beach
            profile_df["id"] = (
                profile_df.beach + "_" + profile_df.section + "_" + profile_df.profile
            )

            # Filter to drop pre-1987 and deep water samples, add to
            # list if profile crosses 0
            profile_df = profile_df[profile_df.date > "1987"]
            profile_df = profile_df[profile_df.z > -3.0]
            if (profile_df.z.min() < datum) & (profile_df.z.max() > datum):
                site_profiles.append(profile_df)

        # If list of profiles contain valid data
        if len(site_profiles) > 0:

            # Combine individual profiles into a single dataframe
            profiles_df = pd.concat(site_profiles)

            # Reproject coords to Albers
            trans = Transformer.from_crs("EPSG:32756", "EPSG:3577", always_xy=True)
            profiles_df["x"], profiles_df["y"] = trans.transform(
                profiles_df.x.values, profiles_df.y.values
            )

            # Compute origin and end points for each profile
            start_xy = profiles_df.groupby(["id"], as_index=False).first()[
                ["id", "x", "y"]
            ]
            end_xy = profiles_df.groupby(["id"], as_index=False).last()[
                ["id", "x", "y"]
            ]
            start_xy = start_xy.rename({"x": "start_x", "y": "start_y"}, axis=1)
            end_xy = end_xy.rename({"x": "end_x", "y": "end_y"}, axis=1)

            # Join origin and end points into dataframe
            profiles_df = pd.merge(left=profiles_df, right=start_xy)
            profiles_df = pd.merge(left=profiles_df, right=end_xy)

            # Compute chainage
            profiles_df["distance"] = profiles_df.apply(
                lambda x: Point(x.start_x, x.start_y).distance(Point(x.x, x.y)), axis=1
            )

            # Drop profiles that have been assigned incorrect profile IDs.
            # To do this, we use a correlation test to determine whether x
            # and y coordinates within each individual profiles fall along a
            # straight line. If a profile has a low correlation (e.g. less
            # than 99.9), it is likely that multiple profile lines have been
            # incorrectly labelled with a single profile ID.
            valid_profiles = lambda x: x[["x", "y"]].corr().abs().iloc[0, 1] > 0.99
            drop = (~profiles_df.groupby("id").apply(valid_profiles)).sum()
            profiles_df = profiles_df.groupby("id").filter(valid_profiles)
            if drop.sum() > 0:
                print(f"\nDropping invalid profiles: {drop:<80}")

            # Restrict profiles to data that falls ocean-ward of the top of
            # the foredune (the highest point in the profile) to remove
            # spurious validation points, e.g. due to a non-shoreline lagoon
            # at the back of the profile
            foredune_dist = (
                profiles_df.groupby(["id", "date"])
                .apply(lambda x: x.distance.loc[x.z.idxmax()])
                .reset_index(name="foredune_dist")
            )
            profiles_df = pd.merge(left=profiles_df, right=foredune_dist)
            profiles_df = profiles_df.loc[
                (profiles_df.distance >= profiles_df.foredune_dist)
            ]

            # Find location and distance to water for datum height (e.g. 0 m AHD)
            intercept_df = (
                profiles_df.groupby(["id", "date"])
                .apply(waterline_intercept, z_val=datum)
                .dropna()
            )

            # If the output contains data
            if len(intercept_df.index) > 0:

                # Join into dataframe
                shoreline_dist = intercept_df.join(
                    profiles_df.groupby(["id", "date"]).first()
                )

                # Compute validation slope and join into dataframe
                slope = val_slope(profiles_df, intercept_df, datum=datum)
                shoreline_dist = shoreline_dist.join(slope.rename("slope"))

                # Keep required columns
                shoreline_dist = shoreline_dist[
                    [
                        "beach",
                        "section",
                        "profile",
                        "name",
                        "source",
                        "foredune_dist",
                        "slope",
                        "start_x",
                        "start_y",
                        "end_x",
                        "end_y",
                        f"{datum}_dist",
                        f"{datum}_x",
                        f"{datum}_y",
                    ]
                ]

                # Export to file
                shoreline_dist.to_csv(fname_out)

    else:
        print(f"Skipping {fname_out:<80}", end="\r")


def preprocess_sadew(fname, datum=0, overwrite=False):

    # Get output filename
    name = Path(fname).stem.split("_")[-1].lower().replace(" ", "")
    fname_out = f"output_data/sadew_{name}.csv"
    print(f"Processing {fname_out:<80}", end="\r")

    # Test if file exists
    if not os.path.exists(fname_out) or overwrite:

        # Load data and set nodata values to NaN
        wide_df = pd.read_csv(fname).replace(-9999, np.nan)
        wide_df = wide_df.loc[:, ~wide_df.columns.str.contains("\.")]
        wide_df.columns = wide_df.columns.str.lower()
        wide_df = wide_df.rename({"easting": "x", "northing": "y"}, axis=1)

        # Determine coordinate transform to use
        utm_54_profiles = [
            730007,
            715055,
            200001,
            710019,
            615010,
            725004,
            725038,
            200047,
            200074,
            200020,
            200031,
            710030,
            620005,
            200046,
            710001,
            620016,
            735007,
            200004,
            710018,
            625002,
            200038,
            730003,
            710006,
            725072,
            715009,
            200128,
            725033,
            200017,
            200126,
            200034,
            200060,
            200030,
            615006,
            710031,
            200053,
            735002,
            200056,
            200049,
            200028,
            200057,
            615003,
            715003,
            735006,
            200127,
            725029,
            725010,
            200025,
            200033,
            200042,
            730010,
            200043,
            200040,
            730004,
            200012,
            200051,
            710023,
            620008,
            725071,
            625003,
            730009,
            615005,
            615007,
            615002,
            200055,
            730006,
            735004,
            200007,
            715056,
            200059,
            625001,
            200008,
            735005,
            715004,
            200054,
            730001,
            725009,
            710022,
            725028,
            730008,
            200019,
            200044,
            200050,
            200032,
            200036,
            710029,
            730002,
            200010,
            615011,
            200052,
            200026,
            710025,
            200021,
            200068,
            200002,
            715006,
            725001,
            200037,
            200013,
            710026,
            620006,
            710027,
            200048,
            620010,
            730005,
            710024,
            200035,
            200006,
            620007,
            710032,
            200122,
            725006,
            201057,
            200005,
            725005,
            710021,
            200129,
            615009,
            715062,
            710017,
            200024,
            735003,
            200045,
            200029,
            620009,
            200039,
            200015,
            200058,
            200124,
            620004,
            620002,
            200041,
            620012,
            625004,
            615004,
            725031,
            200003,
            725008,
            200011,
        ]
        crs = "EPSG:28354" if wide_df.profile[0] in utm_54_profiles else "EPSG:28353"

        # Reproject coords to Albers and create geodataframe
        trans = Transformer.from_crs(crs, "EPSG:3577", always_xy=True)
        wide_df["x"], wide_df["y"] = trans.transform(wide_df.x.values, wide_df.y.values)

        # Reshape into long format with each observation on a new row
        profile_df = pd.melt(
            wide_df.drop("sample_no", axis=1),
            id_vars=["x", "y", "profile"],
            value_name="z",
        ).dropna()

        # Extract date info
        profile_df["date"] = profile_df["variable"].str[1:].str.strip()
        profile_df["date"] = pd.to_datetime(
            profile_df["date"], format="%d%m%Y", errors="coerce"
        )
        profile_df = profile_df.drop("variable", axis=1)

        # Restrict to post 1987 and pre 2020
        profile_df = profile_df[
            (profile_df.date.dt.year > 1987) & (profile_df.date.dt.year < 2020)
        ]

        # Add unique ID column
        profile_df["beach"] = "sadew"
        profile_df["section"] = "all"
        profile_df["profile"] = profile_df["profile"].astype(str)
        profile_df["id"] = (
            profile_df.beach + "_" + profile_df.section + "_" + profile_df.profile
        )
        profile_df["source"] = "hydrographic survey"
        profile_df["name"] = "sadew"

        # Compute origin points for each profile
        profile_df = profile_df.assign(
            start_x=wide_df.iloc[0, 2],
            start_y=wide_df.iloc[0, 3],
            end_x=wide_df.iloc[-1, 2],
            end_y=wide_df.iloc[-1, 3],
        )

        # Compute chainage
        profile_df["distance"] = profile_df.apply(
            lambda x: math.hypot(x.x - x.start_x, x.y - x.start_y), axis=1
        )

        # Find location and distance to water for datum height (e.g. 0 m AHD)
        intercept_df = (
            profile_df.groupby(["id", "date"])
            .apply(waterline_intercept, z_val=datum)
            .dropna()
        )

        # If the output contains data
        if len(intercept_df.index) > 0:

            # Join into dataframe
            shoreline_dist = intercept_df.join(
                profile_df.groupby(["id", "date"]).first()
            )

            # Compute validation slope and join into dataframe
            slope = val_slope(profile_df, intercept_df, datum=datum)
            shoreline_dist = shoreline_dist.join(slope.rename("slope"))

            # Keep required columns
            shoreline_dist = shoreline_dist[
                [
                    "beach",
                    "section",
                    "profile",
                    "name",
                    "source",
                    "slope",
                    "start_x",
                    "start_y",
                    "end_x",
                    "end_y",
                    f"{datum}_dist",
                    f"{datum}_x",
                    f"{datum}_y",
                ]
            ]

            # Export to file
            shoreline_dist.to_csv(fname_out)

    else:
        print(f"Skipping {fname:<80}", end="\r")


def preprocess_sunshinecoast(site, datum=0, overwrite=False):

    # Standardise beach name from site name
    beach = site[2:].replace(" ", "").lower()
    fname_out = f"output_data/sunshinecoast_{beach}.csv"
    print(f"Processing {fname_out:<80}", end="\r")

    # Test if file exists
    if not os.path.exists(fname_out) or overwrite:

        # Obtain list of files
        file_list = glob.glob(
            f"input_data/sunshinecoast/Survey Database correct data/*{site}/**/*.xlsx",
            recursive=True,
        )

        # Output list to hold data
        site_profiles = []

        for i, survey_fname in enumerate(file_list):

            # Load data
            survey_data = pd.read_excel(survey_fname)
            profile_df = survey_data.iloc[3:, :2].astype("float32")
            profile_df.columns = ["distance", "z"]

            # Get data from string
            date = (
                os.path.basename(survey_fname)
                .replace("_", "/")
                .replace("a", "")
                .replace("b", "")
                .replace("00", "01")[0:10]
            )
            profile_df["date"] = pd.to_datetime(date)
            profile_df["beach"] = beach
            profile_df["section"] = "na"
            profile_df["profile"] = survey_fname.split("/")[4]
            profile_df["name"] = "sunshinecoast"
            profile_df["source"] = "hydrographic survey"
            profile_df["id"] = (
                profile_df.beach + "_" + profile_df.section + "_" + profile_df.profile
            )

            # Assign header metadata
            profile_df["start_x"] = survey_data.iloc[2, 2]
            profile_df["start_y"] = survey_data.iloc[2, 3]
            profile_df["bearing"] = survey_data.iloc[1, 3]

            # Fix Kings Beach
            if "KB" in survey_fname.split("/")[4]:
                profile_df["profile"] = survey_fname.split("/")[5]
                profile_df["bearing"] = 125.7
                profile_df["id"] = (
                    profile_df.beach
                    + "_"
                    + profile_df.section
                    + "_"
                    + profile_df.profile
                )

            # Compute
            profile_df["end_y"], profile_df["end_x"] = dist_angle(
                profile_df["start_x"].iloc[0],
                profile_df["start_y"].iloc[0],
                8000,
                profile_df["bearing"].iloc[0],
            )

            # Filter to drop pre-1987 and deep water samples, add to list if any
            # data is available above 0 MSL
            if (
                (profile_df.date.min().year > 1987)
                & (profile_df.z.min() < datum)
                & (profile_df.z.max() > datum)
            ):
                site_profiles.append(profile_df)

        # If list of profiles contain valid data
        if len(site_profiles) > 0:

            # Combine into a single dataframe
            profiles_df = pd.concat(site_profiles)

            # Add coordinates at every supplied distance along transects
            profiles_df[["x", "y"]] = profiles_df.apply(
                lambda x: pd.Series(
                    dist_along_transect(
                        x.distance, x.start_x, x.start_y, x.end_x, x.end_y
                    )
                ),
                axis=1,
            )

            # Convert coordinates to Australian Albers
            trans = Transformer.from_crs("EPSG:32756", "EPSG:3577", always_xy=True)
            profiles_df["start_x"], profiles_df["start_y"] = trans.transform(
                profiles_df["start_x"].values, profiles_df["start_y"].values
            )
            profiles_df["end_x"], profiles_df["end_y"] = trans.transform(
                profiles_df["end_x"].values, profiles_df["end_y"].values
            )
            profiles_df["x"], profiles_df["y"] = trans.transform(
                profiles_df["x"].values, profiles_df["y"].values
            )

            # Readjust distance measurements to Australian Albers. This ensures they
            # are a valid comparison against the Albers-based DEA Coastlines distances
            profiles_df["distance"] = profiles_df.apply(
                lambda x: Point(x.start_x, x.start_y).distance(Point(x.x, x.y)), axis=1
            )

            # Find location and distance to water for datum height (e.g. 0 m AHD)
            intercept_df = (
                profiles_df.groupby(["id", "date"])
                .apply(waterline_intercept, z_val=datum)
                .dropna()
            )

            # If the output contains data
            if len(intercept_df.index) > 0:

                # Join into dataframe
                shoreline_dist = intercept_df.join(
                    profiles_df.groupby(["id", "date"]).agg(
                        lambda x: pd.Series.mode(x).iloc[0]
                    )
                )

                # Compute validation slope and join into dataframe
                slope = val_slope(profiles_df, intercept_df, datum=datum)
                shoreline_dist = shoreline_dist.join(slope.rename("slope"))

                # Keep required columns
                shoreline_dist = shoreline_dist[
                    [
                        "beach",
                        "section",
                        "profile",
                        "name",
                        "source",
                        "slope",
                        "start_x",
                        "start_y",
                        "end_x",
                        "end_y",
                        f"{datum}_dist",
                        f"{datum}_x",
                        f"{datum}_y",
                    ]
                ]

                # Export to file
                shoreline_dist.to_csv(fname_out)

    else:
        print(f"Skipping {fname_out:<80}", end="\r")


def preprocess_tasmarc(site, datum=0, overwrite=True):
    def _tasmarc_metadata(profile):

        # Open file
        with open(profile, "r") as profile_data:

            # Load header data (first 20 rows starting with "#")
            header = takewhile(lambda x: x.startswith(("#", "&", " ")), profile_data)
            header = list(header)[0:19]

            # List of metadata to extract
            meta_list = [
                "LONGITUDE",
                "LATITUDE",
                "START DATE/TIME",
                "TRUE BEARING TRANSECT DEGREES",
                "SURVEY METHOD",
            ]

            # Extract metadata for each metadata type in list
            meta_extract = map(
                lambda m: [
                    row.replace(f"# {m} ", "").strip(" \n")
                    for row in header
                    if m in row
                ][0],
                meta_list,
            )
            lon, lat, date, bearing, source = meta_extract

            # Return metadata
            return {
                "lon": lon,
                "lat": lat,
                "date": date[0:10],
                "bearing": bearing,
                "source": source,
            }

    # List of invalid profiles
    invalid_list = ["shelly_beach_north"]  # incorrect starting coord

    # Set output name
    fname_out = f"output_data/tasmarc_{site}.csv"
    print(f"Processing {fname_out:<80}", end="\r")

    # Test if file exists
    if not os.path.exists(fname_out) or overwrite:

        # List of profile datasets to iterate through
        profile_list = glob.glob(f"input_data/tasmarc/{site}/*.txt")

        # Remove invalid profiles
        profile_list = [
            profile
            for profile in profile_list
            if not any(invalid in profile for invalid in invalid_list)
        ]

        # Output list to hold data
        site_profiles = []

        for profile in profile_list:

            # Load data and remove invalid data
            profile_df = pd.read_csv(
                profile,
                comment="#",
                delim_whitespace=True,
                header=None,
                usecols=[0, 1, 2],
                engine="python",
                names=["distance", "z", "flag"],
            )

            # Remove invalid data
            profile_df = profile_df[profile_df.flag == 2]
            profile_df = profile_df.drop("flag", axis=1)

            # Add metadata from file and coerce to floats
            profile_df = profile_df.assign(**_tasmarc_metadata(profile))
            profile_df = profile_df.apply(pd.to_numeric, errors="ignore")

            # Set to datetime
            profile_df["date"] = pd.to_datetime(profile_df["date"])
            profile_df = profile_df[profile_df.date > "1987"]

            # Set remaining location metadata and ID
            profile_i = os.path.basename(profile).replace(site, "")[1:-15]
            profile_df["profile"] = (
                profile_i.replace("_", "") if len(profile_i) > 0 else "middle"
            )
            profile_df["beach"] = site.replace("_", "")
            profile_df["section"] = "all"
            profile_df["name"] = "tasmarc"
            profile_df["id"] = (
                profile_df.beach + "_" + profile_df.section + "_" + profile_df.profile
            )

            # Remove strings from distance column
            profile_df["distance"] = profile_df.distance.apply(
                pd.to_numeric, errors="coerce"
            )

            # Filter to drop pre-1987 and deep water samples, add to list if any
            # data is available above 0 MSL
            if (profile_df.z.min() < datum) & (profile_df.z.max() > datum):
                site_profiles.append(profile_df)

        # If list of profiles contain valid data
        if len(site_profiles) > 0:

            # Combine into a single dataframe
            profiles_df = pd.concat(site_profiles)

            # Extend survey lines out from start coordinates using supplied angle
            coords_end = profiles_df.apply(
                lambda x: dist_angle(x.lon, x.lat, 0.002, x.bearing), axis=1
            )
            profiles_df = pd.concat([profiles_df, coords_end], axis=1).drop(
                "bearing", axis=1
            )

            # Rename fields
            profiles_df = profiles_df.rename(
                {"lat": "start_y", "lon": "start_x"}, axis=1
            )

            # Reproject coords to Albers and create geodataframe
            trans = Transformer.from_crs("EPSG:4326", "EPSG:3577", always_xy=True)
            profiles_df["start_x"], profiles_df["start_y"] = trans.transform(
                profiles_df.start_x.values, profiles_df.start_y.values
            )
            profiles_df["end_x"], profiles_df["end_y"] = trans.transform(
                profiles_df.end_x.values, profiles_df.end_y.values
            )

            # Add coordinates for every distance along transects
            profiles_df[["x", "y"]] = profiles_df.apply(
                lambda x: pd.Series(
                    dist_along_transect(
                        x.distance, x.start_x, x.start_y, x.end_x, x.end_y
                    )
                ),
                axis=1,
            )

            # Find location and distance to water for datum height (0 m AHD)
            intercept_df = (
                profiles_df.groupby(["id", "date"])
                .apply(waterline_intercept, z_val=datum)
                .dropna()
            )

            # If the output contains data
            if len(intercept_df.index) > 0:

                # Join into dataframe
                shoreline_dist = intercept_df.join(
                    profiles_df.groupby(["id", "date"]).agg(
                        lambda x: pd.Series.mode(x).iloc[0]
                    )
                )

                # Compute validation slope and join into dataframe
                slope = val_slope(profiles_df, intercept_df, datum=datum)
                shoreline_dist = shoreline_dist.join(slope.rename("slope"))

                # Keep required columns
                shoreline_dist = shoreline_dist[
                    [
                        "beach",
                        "section",
                        "profile",
                        "name",
                        "source",
                        "slope",
                        "start_x",
                        "start_y",
                        "end_x",
                        "end_y",
                        f"{datum}_dist",
                        f"{datum}_x",
                        f"{datum}_y",
                    ]
                ]

                # Standardise source column
                standardise_source(shoreline_dist)

                # Export to file
                shoreline_dist.to_csv(fname_out)

        else:
            print(f"Skipping {fname_out:<80}", end="\r")


def preprocess_moruya(fname_out, datum=0, overwrite=False):
    def azimuth(point1, point2):
        """azimuth between 2 shapely points (interval 0 - 360)"""
        angle = np.arctan2(point2.x - point1.x, point2.y - point1.y)
        return np.degrees(angle) if angle > 0 else np.degrees(angle) + 360

    # Test if file exists
    if not os.path.exists(fname_out) or overwrite:

        # Import data
        profiles_wide = pd.read_csv("input_data/moruya/moruya_processed.csv")

        # Set to datetime
        profiles_wide["date"] = pd.to_datetime(profiles_wide["date"], format="%d/%m/%Y")

        # Reshape to long format with one row per measurement
        profiles_df = profiles_wide.melt(
            id_vars=profiles_wide.columns[0:8],
            value_vars=profiles_wide.columns[8:],
            var_name="distance",
            value_name="z",
        )

        # Drop rows with no elevation measurements
        profiles_df = profiles_df.dropna(axis=0, subset=["z"])

        # Convert units to metres
        profiles_df["z"] = profiles_df["z"] * 0.01

        # Convert distance strings to numeric distances
        profiles_df["distance"] = profiles_df.apply(
            lambda x: int(x.distance[5:]), axis=1
        )

        # Generate survey lines start and end points
        profiles_df["angle"] = profiles_df.apply(
            lambda x: azimuth(Point(x.start_x, x.start_y), Point(x.mid_x, x.mid_y)),
            axis=1,
        )
        profiles_df[["end_y", "end_x"]] = profiles_df.apply(
            lambda x: dist_angle(x.start_x, x.start_y, 1000, x.angle), axis=1
        )
        profiles_df = profiles_df.drop(["mid_y", "mid_x"], axis=1)

        # Reproject coords to Albers and create geodataframe
        trans = Transformer.from_crs("EPSG:32756", "EPSG:3577", always_xy=True)
        profiles_df["start_x"], profiles_df["start_y"] = trans.transform(
            profiles_df.start_x.values, profiles_df.start_y.values
        )
        profiles_df["end_x"], profiles_df["end_y"] = trans.transform(
            profiles_df.end_x.values, profiles_df.end_y.values
        )

        # Create site and source column
        profiles_df["id"] = (
            profiles_df.beach
            + "_"
            + profiles_df.section
            + "_"
            + profiles_df.profile.map(str)
        )
        profiles_df["source"] = "emery/levelling"
        profiles_df["name"] = "moruya"

        # Add coordinates at every supplied distance along transects
        profiles_df[["x", "y"]] = profiles_df.apply(
            lambda x: pd.Series(
                dist_along_transect(x.distance, x.start_x, x.start_y, x.end_x, x.end_y)
            ),
            axis=1,
        )

        # Find location and distance to water for datum height (e.g. 0 m AHD)
        intercept_df = (
            profiles_df.groupby(["id", "date"])
            .apply(waterline_intercept, z_val=datum)
            .dropna()
        )

        # If the output contains data
        if len(intercept_df.index) > 0:

            # Join into dataframe
            shoreline_dist = intercept_df.join(
                profiles_df.groupby(["id", "date"]).agg(
                    lambda x: pd.Series.mode(x).iloc[0]
                )
            )

            # Compute validation slope and join into dataframe
            slope = val_slope(profiles_df, intercept_df, datum=datum)
            shoreline_dist = shoreline_dist.join(slope.rename("slope"))

            # Keep required columns
            shoreline_dist = shoreline_dist[
                [
                    "beach",
                    "section",
                    "profile",
                    "name",
                    "source",
                    "slope",
                    "start_x",
                    "start_y",
                    "end_x",
                    "end_y",
                    f"{datum}_dist",
                    f"{datum}_x",
                    f"{datum}_y",
                ]
            ]

            # Remove dodgy Pedro 4 transect
            shoreline_dist = shoreline_dist.drop("pedro_all_4")

            # Export to file
            shoreline_dist.to_csv(fname_out)


def waterbody_mask(input_data, modification_data, bbox):
    """
    Generates a raster mask for DEACoastLines based on the
    SurfaceHydrologyPolygonsRegional.gdb dataset, and a vector
    file containing minor modifications to this dataset (e.g.
    features to remove or add to the dataset).

    The mask returns True for perennial 'Lake' features, any
    'Aquaculture Area', 'Estuary', 'Watercourse Area', 'Salt
    Evaporator', and 'Settling Pond' features. Features of
    type 'add' from the modification data file are added to the
    mask, while features of type 'remove' are removed.
    """

    # Import SurfaceHydrologyPolygonsRegional data
    waterbody_gdf = gpd.read_file(input_data, bbox=bbox).to_crs(bbox.crs)

    # Restrict to coastal features
    lakes_bool = (waterbody_gdf.FEATURETYPE == "Lake") & (
        waterbody_gdf.PERENNIALITY == "Perennial"
    )
    other_bool = waterbody_gdf.FEATURETYPE.isin(
        [
            "Aquaculture Area",
            "Estuary",
            "Watercourse Area",
            "Salt Evaporator",
            "Settling Pond",
        ]
    )
    waterbody_gdf = waterbody_gdf[lakes_bool | other_bool]

    # Load in modification dataset and select features to remove/add
    mod_gdf = gpd.read_file(modification_data, bbox=bbox).to_crs(bbox.crs)
    to_remove = mod_gdf[mod_gdf["type"] == "remove"]
    to_add = mod_gdf[mod_gdf["type"] == "add"]

    # Remove and add features
    if len(to_remove.index) > 0:
        if len(waterbody_gdf.index) > 0:
            waterbody_gdf = gpd.overlay(waterbody_gdf, to_remove, how="difference")
    if len(to_add.index) > 0:
        if len(waterbody_gdf.index) > 0:
            waterbody_gdf = gpd.overlay(waterbody_gdf, to_add, how="union")
        else:
            waterbody_gdf = to_add

    return waterbody_gdf


def smartline_attrs(val_gdf, bbox, advanced=True):
    def nearest_features(input_gdf, comp_gdf, cols="erode_v"):

        # Verify that both datasets are in the same CRS
        assert input_gdf.crs == comp_gdf.crs

        i_mins = []
        dist_mins = []

        # Loop through each input and comparison feature, and create lists of
        # distances and nearest feature indices
        for input_ft in input_gdf.geometry:
            dists = [input_ft.distance(comp_ft) for comp_ft in comp_gdf.geometry]
            dist_mins.append(np.min(dists))
            i_mins.append(np.argmin(dists))

        # Use indices to extract relevant rows from comparison data,
        # and add a column giving distance to these nearest features
        out_gdf = comp_gdf.iloc[i_mins][[cols] if type(cols) == str else cols]
        out_gdf["near_dist"] = dist_mins

        # Verify that there is one output for every input feature
        assert len(input_gdf) == input_gdf.shape[0]

        return out_gdf

    # Load Smartline advanced data
    smartline = (
        gpd.read_file(
            "https://dea-public-data.s3.ap-southeast-2.amazonaws.com/derivative/dea_coastlines/supplementary/Smartline.gpkg",
            bbox=bbox.buffer(100),
        )
        .to_crs("EPSG:3577")
        .rename({"INTERTD1_V": "smartline"}, axis=1)
    )

    # Identify unique profiles (no need to repeat Smartline extraction for each time)
    inds = [i[1][0] for i in val_gdf.groupby("id").groups.items()]
    input_gdf = val_gdf.iloc[inds].reset_index()

    # Find nearest features
    nearest_df = nearest_features(
        input_gdf=input_gdf, comp_gdf=smartline, cols="smartline"
    ).reset_index()

    # # Spatial join unique profiles to smartline and drop unclassified
    input_gdf["smartline"] = nearest_df.smartline
    smartline_df = input_gdf[["id", "smartline"]]
    smartline_df = smartline_df.loc[smartline_df.smartline != "Unclassified"]

    return smartline_df


def deacl_validation(
    val_path,
    deacl_path,
    datum=0,
    prefix="temp",
    overwrite=False,
    sat_label="deacl",
    val_label="val",
    layer_name="shorelines_annual",
):

    # Set up output file name
    out_name = val_path.split("/")[-1]
    out_name = f"data/validation/processed/outputs_{prefix}_{out_name}"

    # Run analysis if file doesn't exist and overwrite is True
    if not os.path.exists(out_name) or overwrite:

        # Load validation data and set section/beach
        print(f"{val_path:<80}", end="\r")
        val_df = pd.read_csv(val_path, parse_dates=["date"])
        cat_cols = ["beach", "section", "profile", "source", "name"]
        val_df[cat_cols] = val_df[cat_cols].astype("str")

        # Get bounding box to load data for
        minx, maxy = val_df.min().loc[[f"{datum}_x", f"{datum}_y"]]
        maxx, miny = val_df.max().loc[[f"{datum}_x", f"{datum}_y"]]
        bbox = gpd.GeoSeries(box(minx, miny, maxx, maxy), crs="EPSG:3577")

        # Read in DEA Coastlines data
        deacl_gdf = (
            gpd.read_file(deacl_path, bbox=bbox.buffer(100), layer=layer_name)
            .to_crs("EPSG:3577")
            .dissolve("year")
            .reset_index()
        )

        # This will fail if no DEA Coastlines data exists for the area
        if (len(deacl_gdf.index) > 0) & (len(val_df.index) > 0):

            # Add year columns and set dtype to allow merging
            deacl_gdf["year"] = deacl_gdf.year.astype("int64")
            val_df["year"] = val_df.date.dt.year

            # Aggregate by year for each ID and compute modal, count and
            # median stats
            modal_vals = val_df.groupby(["year", "id"])[cat_cols].agg(
                lambda x: pd.Series.mode(x).iloc[0]
            )
            count_vals = val_df.groupby(["year", "id"]).year.count().rename("n")
            median_vals = val_df.groupby(["year", "id"]).median()

            # Combine all aggregated stats into one dataframe
            val_df = pd.concat(
                [median_vals, modal_vals, count_vals], axis=1
            ).reset_index()

            # Convert validation start and end locations to linestrings
            val_geometry = val_df.apply(
                lambda x: LineString([(x.start_x, x.start_y), (x.end_x, x.end_y)]),
                axis=1,
            )

            # Convert geometries to GeoDataFrame
            val_gdf = gpd.GeoDataFrame(
                data=val_df, geometry=val_geometry, crs="EPSG:3577"
            ).reset_index()

            # Join Smartline data
            val_gdf = val_gdf.merge(
                smartline_attrs(val_gdf, bbox.buffer(1000)), how="left", on="id"
            )

            # Match each shoreline contour to each date in validation data
            results_df = val_gdf.merge(
                deacl_gdf, on="year", suffixes=("_val", "_deacl")
            )

            # For each row, identify where profile intersects with waterline
            results_df["intersect"] = results_df.apply(
                lambda x: x.geometry_val.intersection(x.geometry_deacl), axis=1
            )

            # Drop any multipart geometries as these are invalid comparisons
            results_df = results_df[
                results_df.apply(lambda x: x.intersect.type == "Point", axis=1)
            ]
            results_df[f"{sat_label}_x"] = gpd.GeoSeries(results_df["intersect"]).x
            results_df[f"{sat_label}_y"] = gpd.GeoSeries(results_df["intersect"]).y

            # For each row, compute distance between origin and intersect
            results_df[f"{sat_label}_dist"] = results_df.apply(
                lambda x: x.intersect.distance(Point(x.start_x, x.start_y)), axis=1
            )

            # If data contains a foredune distance field, drop invalid
            # validation points where DEA CoastLines intersection occurs
            # behind the foredune
            if "foredune_dist" in results_df:
                valid = results_df[f"{sat_label}_dist"] >= results_df.foredune_dist
                results_df = results_df.loc[valid]

            # TO REMOVE
            if "slope" not in results_df:
                results_df["slope"] = np.nan

            # If enough data is returned:
            if len(results_df.index) > 0:

                # Rename for consistency
                results_df = results_df.rename(
                    {
                        f"{datum}_dist": f"{val_label}_dist",
                        f"{datum}_x": f"{val_label}_x",
                        f"{datum}_y": f"{val_label}_y",
                    },
                    axis=1,
                )

                # Calculate difference
                results_df["error_m"] = results_df.val_dist - results_df.deacl_dist

                # Add lon and lat columns
                centroid_points = gpd.GeoSeries(
                    results_df.geometry_val
                ).centroid.to_crs("EPSG:4326")
                results_df["lon"] = centroid_points.x.round(3)
                results_df["lat"] = centroid_points.y.round(3)

                # Keep correct columns
                results_df = results_df[
                    [
                        "id",
                        "year",
                        "beach",
                        "section",
                        "profile",
                        "name",
                        "source",
                        "certainty",
                        "n",
                        "lon",
                        "lat",
                        "slope",
                        "smartline",
                        "start_x",
                        "start_y",
                        "end_x",
                        "end_y",
                        f"{val_label}_x",
                        f"{val_label}_y",
                        f"{val_label}_dist",
                        f"{sat_label}_x",
                        f"{sat_label}_y",
                        f"{sat_label}_dist",
                        "error_m",
                    ]
                ]

                # Export data
                results_df.to_csv(out_name, index=False)

    # Else skip
    else:
        print("Skipping; file either exists or overwrite set to false")


@click.command()
@click.option(
    "--inputs_path",
    type=str,
    required=True,
    help="",
)
@click.option(
    "--deacl_path",
    type=str,
    required=True,
    help="",
)
@click.option(
    "--prefix",
    type=str,
    default="temp",
    help="",
)
@click.option(
    "--datum",
    type=str,
    default=0,
    help="",
)
@click.option(
    "--overwrite",
    type=bool,
    default=False,
    help="",
)
@click.option(
    "--layer_name",
    type=str,
    default="shorelines_annual",
    help="",
)
@click.option(
    "--append_stats",
    type=bool,
    default=False,
    help="",
)
@click.option(
    "--parallelised",
    type=bool,
    default=True,
    help="",
)
@click.option(
    "--markdown_report",
    type=bool,
    default=False,
    help="",
)
def validation_cli(
    inputs_path,
    deacl_path,
    prefix,
    datum,
    overwrite,
    layer_name,
    append_stats,
    parallelised,
    markdown_report,
):

    # Simplify smartline categories
    rename_dict = {
        "Beachrock undiff": "rocky",
        "Beachrock undiff dominant": "rocky",
        "Boulder or shingle-grade beach undiff": "rocky",
        "Boulder groyne or breakwater undiff": "rocky",
        "Flat boulder deposit (rock) undiff": "rocky",
        "Hard bedrock shore": "rocky",
        "Hard bedrock shore inferred": "rocky",
        "Hard rock cliff (>5m)": "rocky",
        "Hard rocky shore platform": "rocky",
        "Rocky shore platform (undiff)": "rocky",
        "Sloping boulder deposit (rock) undiff": "rocky",
        "Sloping hard rock shore": "rocky",
        "Sloping soft `bedrock¿ shore": "rocky",
        "Sloping soft \u2018bedrock\u2019 shore": "rocky",
        "Soft `bedrock¿ shore inferred": "rocky",
        "Soft `bedrock¿ shore platform": "rocky",
        "Beach (sediment type undiff)": "sandy",
        "Fine-medium sand beach": "sandy",
        "Fine-medium sandy tidal flats": "sandy",
        "Mixed sand and shell beach": "sandy",
        "Mixed sandy shore undiff": "sandy",
        "Perched sandy beach (undiff)": "sandy",
        "Sandy beach undiff": "sandy",
        "Sandy beach with cobbles/pebbles (rock)": "sandy",
        "Sandy shore undiff": "sandy",
        "Sandy tidal flats": "sandy",
        "Sandy tidal flats with coarse stony debris": "sandy",
        "Sandy tidal flats, no bedrock protruding": "sandy",
        "Sloping coffee rock deposit": "rocky",
        "Muddy tidal flats": "muddy",
        "Tidal flats (sediment undiff)": "muddy",
        "Artificial shoreline undiff": "rocky",
        "Artificial boulder structures undiff": "rocky",
        "Boulder revetment": "rocky",
        "Boulder seawall": "rocky",
        "Concrete sea wall": "rocky",
        "Piles (Jetty)": "rocky",
        "Coarse sand beach": "sandy",
    }

    # Find data to process
    val_paths = glob.glob(f"{inputs_path}*.csv")

    if parallelised:

        from concurrent.futures import ProcessPoolExecutor
        from tqdm import tqdm
        from itertools import repeat

        args = [deacl_path, datum, prefix, overwrite, "deacl", "val", layer_name]

        with ProcessPoolExecutor() as executor:

            # Apply func in parallel
            groups = val_paths
            to_iterate = (groups, *(repeat(i, len(groups)) for i in args))
            tqdm(executor.map(deacl_validation, *to_iterate), total=len(groups))

    else:

        # Run validation for each input dataset
        for val_path in val_paths:

            deacl_validation(
                val_path=val_path,
                deacl_path=deacl_path,
                datum=datum,
                prefix=prefix,
                overwrite=overwrite,
                layer_name=layer_name,
            )

    print("Combining data")
    outputs_list = glob.glob(f"data/validation/processed/outputs_{prefix}_*.csv")
    outputs_df = pd.concat([pd.read_csv(csv) for csv in outputs_list])

    # Rename smartline categories to smaller subset
    outputs_df["smartline"] = outputs_df.smartline.replace(rename_dict)

    # Run stats
    stats_df = deacl_val_stats(
        val_dist=outputs_df.val_dist,
        deacl_dist=outputs_df.deacl_dist,
        n=outputs_df.n,
        remove_bias=False,
    )

    # Transpose and add index time and prefix name
    stats_df = pd.DataFrame({pd.to_datetime("now"): stats_df}).T.assign(name=prefix)
    stats_df.index.name = "time"
    filename = f"data/validation/processed/stats_{prefix}.csv"

    if append_stats:

        stats_df.to_csv(
            filename,
            mode="a",
            header=(not os.path.exists(filename)),
        )

        # Read in stats from disk (ensures we get older results)
        stats_df = pd.read_csv(filename, index_col=0, parse_dates=True)

    else:

        stats_df.to_csv(filename)

    # Plot data
    fig, axes = plt.subplots(1, 2, figsize=(15, 7.5))

    # Extract integration test run times and convert to local time
    times_local = stats_df.index.tz_localize("UTC").tz_convert(tz="Australia/Canberra")
    stats_df.index = times_local

    # Compute stats
    n, mae, rmse, stdev, corr, bias, _ = stats_df.iloc[-1]
    offset_str = "(towards land)" if bias > 0 else "(towards ocean)"

    # Plot latest integration test results as scatterplot
    outputs_df.plot.scatter(ax=axes[0], x="val_dist", y="deacl_dist")
    axes[0].plot([0, 150], [0, 150], linestyle="--", color="black", linewidth=0.5)
    axes[0].set_title(
        "Latest integration test run:\nDEA Coastlines vs. validation scatterplot"
    )
    axes[0].set_ylabel("DEA Coastlines shoreline positions (m)")
    axes[0].set_xlabel("Validation shoreline positions (m)")
    axes[0].annotate(
        f"Mean Absolute Error: {mae:.2f} m\n"
        f"RMSE: {rmse:.2f} m\n"
        f"Standard deviation: {stdev:.2f} m\n"
        f"Bias: {bias:.2f} m {offset_str}\n"
        f"Correlation: {corr:.3f}\n",
        xy=(0.05, 0.95),
        horizontalalignment="left",
        verticalalignment="top",
        xycoords="axes fraction",
    )

    # Plot all integration test accuracies and biases over time
    stats_df.rmse.rename("RMSE").plot(ax=axes[1], style=".-", legend=True)
    min_q, max_q = stats_df.rmse.quantile((0.1, 0.9)).values
    axes[1].fill_between(stats_df.index, min_q, max_q, alpha=0.2)

    stats_df.mae.rename("MAE").plot(ax=axes[1], style=".-", legend=True)
    min_q, max_q = stats_df.mae.quantile((0.1, 0.9)).values
    axes[1].fill_between(stats_df.index, min_q, max_q, alpha=0.2)

    stats_df.bias.rename("Bias").plot(ax=axes[1], style=".-", legend=True)
    min_q, max_q = stats_df.bias.quantile((0.1, 0.9)).values
    axes[1].fill_between(stats_df.index, min_q, max_q, alpha=0.2)
    axes[1].set_title("Accuracy and bias across\n all integration test runs")
    axes[1].set_ylabel("Metres (m)")
    axes[1].set_xlabel(None)

    # Add overall title
    plt.suptitle(
        f"Latest DEA Coastlines integration test outputs validation ({str(times_local[-1])[0:16]})",
        size=14,
        fontweight="bold",
        y=1.03,
    )

    # Export to file
    plt.savefig(f"data/validation/processed/stats_{prefix}.png", bbox_inches="tight")

    if markdown_report:

        # Create markdown report
        from mdutils.mdutils import MdUtils
        from mdutils import Html

        # Calculate recent change and convert to plain text
        stats_df_temp = stats_df.drop("name", axis=1).copy()
        stats_df_temp['bias'] = stats_df_temp['bias'].abs()
        recent_diff = stats_df_temp.diff(1).iloc[-1].to_frame("diff")
        recent_diff.loc['corr'] = -recent_diff.loc['corr']  # Invert as higher corrs are good
        recent_diff.loc[recent_diff["diff"] < 0, "prefix"] = ":heavy_check_mark: improved by "
        recent_diff.loc[recent_diff["diff"] == 0, "prefix"] = ":heavy_minus_sign: no change"
        recent_diff.loc[recent_diff["diff"] > 0, "prefix"] = ":heavy_exclamation_mark: worsened by "
        recent_diff["suffix"] = recent_diff["diff"].abs().round(3).replace({0: ""})
        recent_diff = recent_diff.prefix.astype(str) + recent_diff.suffix.astype(str).str[0:5]
        
        mdFile = MdUtils(file_name="tests/README.md", title="Integration tests")
        mdFile.new_paragraph(
            "> *This readme is automatically generated by the ``coastlines.validation.py`` script; edits should be made to the ``validation_cli`` function [located here](../coastlines/validation.py).*"
        )
        mdFile.new_paragraph(
            "This directory contains integration tests that are run to verify that the entire Coastlines code runs correctly. The ``test_coastline.py`` file runs a simple Coastlines analysis over a single beach (Narrabeen Beach in northern Sydney), using the DEA Coastlines [Command Line Interface (CLI) tools](../notebooks/DEACoastlines_generation_CLI.ipynb) to run the raster, vector, continental layers and validation analysis steps."
        )

        mdFile.new_header(level=1, title="Validation")
        mdFile.new_paragraph(
            "In addition to testing whether the code runs without errors, we also run a small-scale validation of the results of the integration tests by comparing them to validation data from the [Narrabeen-Collaroy Beach Survey Program](https://doi.org/10.1038/sdata.2016.24). The ensures that the code both works, and generates sensible results."
        )
        mdFile.new_paragraph(
            "> Note that this integration test validation is for a single site and a limited number of years; for a full validation of the DEA Coastlines product, refer to [Bishop-Taylor et al. 2021](https://doi.org/10.1016/j.rse.2021.112734)."
        )

        mdFile.new_header(level=2, title="Latest integration test validation results")
        mdFile.new_paragraph(
            f"The latest integration test completed at **{str(stats_df.index[-1])[0:16]}**. "
            f"Compared to the previous run, it had an:"
        )
        items = [f"RMSE accuracy of **{stats_df.rmse[-1]:.2f} m ({recent_diff.rmse})**", 
                 f"MAE accuracy of **{stats_df.mae[-1]:.2f} m ({recent_diff.mae})**", 
                 f"Bias of **{stats_df['bias'][-1]:.2f} m ({recent_diff['bias']})**",
                 f"Pearson correlation of **{stats_df['corr'][-1]:.3f} ({recent_diff['corr']})**"]
        mdFile.new_list(items=items)
        mdFile.new_paragraph(Html.image(path=f"stats_tests.png", size="950"))
        mdFile.create_md_file()


if __name__ == "__main__":
    validation_cli()
