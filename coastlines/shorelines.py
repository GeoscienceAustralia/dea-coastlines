import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio.features
import xarray as xr
from shapely.geometry import shape

from coastlines.timeseries import _block_theil_sen

def vectorise_uncertainty(
    uncertainty_stack: xr.DataArray,
    count_stack: xr.DataArray,
    uncertainty_threshold: float = 10.0,
    count_threshold: float = 3,
    sieve_size: int = 10
) -> dict[str, gpd.GeoDataFrame]:
    """
    Converts timeseries uncertainty and count rasters into a dictionary of 
    dissolved vector polygon masks indicating 'good', 'high uncertainty', 
    and 'low observations' areas per timestep.
    
    Parameters
    ----------
    uncertainty_stack : xr.DataArray
        A 3D array of spatial uncertainty values with a 'time' dimension.
    count_stack : xr.DataArray
        A 3D array of observation counts with a 'time' dimension.
    uncertainty_threshold : float, optional
        Pixels with uncertainty above or equal to this threshold will be
        flagged as high uncertainty. Default is 10.0.
    count_threshold : float, optional
        Pixels with counts below or equal to this threshold will be
        flagged as low counts. Default is 3.
    sieve_size : int, optional
        The minimum number of connected pixels required to retain an isolated 
        patch of certainty/uncertainty. Default is 10.
        
    Returns
    -------
    dict
        A dictionary mapping each timestep (as a string) to a GeoDataFrame 
        containing dissolved certainty polygons.
    """
    certainty_polygons = {}
    
    # Merge into a dataset to guarantee time dimensions are perfectly aligned during iteration
    ds = xr.Dataset({
        "uncertainty": uncertainty_stack, 
        "count": count_stack
    })
    
    for time_val, group in ds.groupby("time"):
        
        unc_da = group["uncertainty"]
        cnt_da = group["count"]
        
        # Create base mask: good data is zero, high uncertainty is one
        raster_mask = xr.where(unc_da >= uncertainty_threshold, 1, 0)
        
        # Apply count mask: low observations become two, overriding the uncertainty mask
        raster_mask = xr.where(cnt_da <= count_threshold, 2, raster_mask).astype(np.int16)
        
        if sieve_size > 0:
            # rasterio.features.sieve supports multi-class integers seamlessly
            sieved_mask = rasterio.features.sieve(
                raster_mask.values, 
                size=sieve_size
            )
        else:
            sieved_mask = raster_mask.values
            
        shapes = rasterio.features.shapes(
            sieved_mask,
            transform=uncertainty_stack.odc.transform
        )
        
        records = [{"geometry": shape(geom), "certainty": val} for geom, val in shapes]
        gdf = gpd.GeoDataFrame(records, crs=uncertainty_stack.odc.crs)
        
        # Remap integer values to descriptive strings
        gdf["certainty"] = gdf["certainty"].map({
            0: "good", 
            1: "high uncertainty", 
            2: "low observations"
        })
        
        gdf = gdf.dissolve(by="certainty").reset_index()
        
        time_str = pd.to_datetime(time_val).strftime("%Y-%m-%d")
        certainty_polygons[time_str] = gdf
        
    return certainty_polygons


def shoreline_uncertainty(
    shorelines_gdf: gpd.GeoDataFrame, 
    certainty_polygons: dict[str, gpd.GeoDataFrame]
) -> gpd.GeoDataFrame:
    """
    Slices a timeseries of continuous shoreline features using corresponding 
    monthly certainty polygons, assigning certainty classes to each segment.
    
    Parameters
    ----------
    shorelines_gdf : gpd.GeoDataFrame
        The continuous, unmasked shoreline vector features containing a 'time' column.
    certainty_polygons : dict
        A dictionary of vectorized polygon masks keyed by time string.
        
    Returns
    -------
    gpd.GeoDataFrame
        A single concatenated GeoDataFrame where shorelines have been split 
        at polygon boundaries and attributed with a 'certainty' column.
    """
    out_list = []
    
    # Group shorelines by their timestamp
    for time_val, contour_gdf in shorelines_gdf.groupby("time"):
        
        # Match the shoreline timestamp format to the polygon dictionary keys
        # Adjust this slice if your shorelines use a different time format
        time_str = str(time_val)[:10]
        
        if time_str in certainty_polygons:
            # Overlay splits the line geometries at the polygon boundaries 
            attributed_lines = contour_gdf.overlay(
                certainty_polygons[time_str], 
                how="intersection"
            )
            out_list.append(attributed_lines)
        else:
            out_list.append(contour_gdf)
            
    # Combine all months back into a single continuous dataset
    return pd.concat(out_list, ignore_index=True)


def rate_change_theilsen(
    sdf_stack: xr.DataArray,
    alpha: float = 0.90
) -> xr.Dataset:
    """
    Calculates the robust temporal rate of shoreline change.

    Parameters
    ----------
    sdf_stack : xr.DataArray
        A 3D array (time, y, x) of tide-corrected signed distance fields.
    alpha : float, optional
        Confidence degree between 0 and 1. Default is 0.90 for a 90%
        confidence interval.

    Returns
    -------
    xr.Dataset
        A dataset containing 2D maps: median rate of change (m/yr),
        upper confidence bound, lower confidence bound, intercept,
        rate uncertainty and significance. Positive rates indicate
        accretion, negative rates indicate erosion.
    """
    # Drop dates that are completely obscured by cloud or no-data,
    # and invert to match DEA Coastlines change conventions
    sdf_stack = -sdf_stack.dropna(dim="time", how="all")

    # Convert datetime objects to continuous decimal years 
    time_years = sdf_stack.time.dt.decimal_year

    # Run the generic Theil-Sen regression to return temporal slopes
    raw_m, raw_m_l, raw_m_h, intercept = xr.apply_ufunc(
        _block_theil_sen,
        sdf_stack,
        time_years,
        kwargs={'alpha': alpha},
        input_core_dims=[["time"], ["time"]],
        output_core_dims=[[], [], [], []],
        vectorize=False,
        dask="parallelized",
        output_dtypes=[np.float32, np.float32, np.float32, np.float32],
    )

    # Temporal slopes represent rates of coastal change (m/yr)
    ds_out = xr.Dataset(
        {
            "rate_optimal": raw_m,
            "rate_low": raw_m_l,
            "rate_high": raw_m_h,
            "rate_intercept": intercept,
        }
    )

    # Calculate final uncertainty bounds from confidence limits
    ds_out["rate_uncertainty"] = np.abs(ds_out.rate_high - ds_out.rate_low)

    # Add significance; neg * neg AND pos * pos == positive
    ds_out["rate_significance"] = (ds_out.rate_low * ds_out.rate_high) > 0

    return ds_out