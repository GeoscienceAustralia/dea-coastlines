import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import linregress

from coastlines.timeseries import (
    # _block_theil_sen, 
    beach_slope_theilsen,
    _slope_dispersion, 
    _extract_optimal_slopes,
    # flag_boomerang_outliers
)


def _calculate_r2(x_da: xr.DataArray, y_da: xr.DataArray) -> float:
    """
    Calculates the coefficient of determination (R-squared) between two 
    DataArrays, safely dropping missing values.
    """
    valid_mask = x_da.notnull() & y_da.notnull()
    x_valid = x_da.where(valid_mask, drop=True).values
    y_valid = y_da.where(valid_mask, drop=True).values
    
    if len(x_valid) < 2:
        return np.nan
        
    return linregress(x_valid, y_valid).rvalue ** 2


def plot_pixel_exploration(
    sdf_stack: xr.DataArray, 
    tide_heights: xr.DataArray, 
    x: float, 
    y: float,
    mask_outliers: bool = False,
    mad_multiplier: int = 6,
    forward_fill: int | None = None,
    color_by_time: bool = False,
    max_dist : int | None = None,
) -> None:
    """
    Extracts, transforms, and visualizes raw, differenced, and detrended 
    shoreline and tide data for a specific spatial coordinate.

    Parameters
    ----------
    sdf_stack : xr.DataArray
        Array containing signed distance fields with dimensions of time, y, and x.
    tide_heights : xr.DataArray
        Array containing tide heights. Can be 1D (time) or 3D (time, y, x).
    x : float
        The longitude or easting of the target pixel.
    y : float
        The latitude or northing of the target pixel.
    mask_outliers : bool, optional
        If True, applies the boomerang outlier flag before differencing, 
        by default False.
    mad_multiplier : int, optional
        Multiplier for the median absolute deviation used in outlier detection,
        by default 6.
    forward_fill : int, optional
        Limit for forward-filling missing differences, by default None.
    color_by_time : bool, optional
        If True, colors all data points by decimal year to visualize temporal 
        trends across all plots. Default is False.
    max_dist : int, optional
        Maximum SDF distance to treat as valid. Can be useful for excluding outliers.
    """
    # Select specific point
    tide_point = tide_heights.sel(x=x, y=y, method="nearest")
    sdf_point = sdf_stack.sel(x=x, y=y, method="nearest")
    valid_data = sdf_point.notnull()

    # Mask out values greater than absolute max_dist
    if max_dist is not None:
        sdf_point = sdf_point.where(np.abs(sdf_point) < max_dist)

    if forward_fill is not None:
        # Forward fill NaN gaps in both SDF and tide data
        sdf_point = sdf_point.ffill(dim="time", limit=forward_fill)
        tide_point = tide_point.where(valid_data).ffill(dim="time", limit=forward_fill)

        # Compute differences
        diff_sdf = sdf_point.diff(dim="time").where(valid_data)
        diff_tides = tide_point.diff(dim="time").where(valid_data)

        # Re-apply valid data mask to remove forward filled data
        # but retain diffs
        sdf_point = sdf_point.where(valid_data)
        tide_point = tide_point.where(valid_data)
    else:
        # Just compute differences directly
        diff_sdf = sdf_point.diff(dim="time")
        diff_tides = tide_point.diff(dim="time")  

    # if mask_outliers:
    #     outlier_mask = flag_boomerang_outliers(
    #         sdf_stack=sdf_point,
    #         mad_multiplier=mad_multiplier,
    #     )
    #     sdf_point = sdf_point.where(~outlier_mask)

    # sdf_point, tide_point = xr.align(
    #     sdf_point.dropna(dim="time", how="all"),
    #     tide_point.dropna(dim="time", how="all"),
    #     join="inner",
    # )
    # diff_sdf, diff_tides = xr.align(
    #     diff_sdf.dropna(dim="time", how="all"),
    #     diff_tides.dropna(dim="time", how="all"),
    #     join="inner",
    # )

    # Extract NumPy arrays once for cleaner syntax downstream
    x_raw, y_raw = tide_point.values, sdf_point.values
    x_diff, y_diff = diff_tides.values, diff_sdf.values
    t_raw = sdf_point.time.values
    t_diff = diff_sdf.time.values

    r2_raw = _calculate_r2(tide_point, sdf_point)
    r2_diff = _calculate_r2(diff_tides, diff_sdf)

    # Calculate Theil-Sen regression for the scatter trendlines
    ts_raw = stats.theilslopes(y_raw, x_raw)
    ts_diff = stats.theilslopes(y_diff, x_diff)

    # Convert regression slopes to physical beach slopes (tan beta)
    raw_slope = 1.0 / ts_raw[0]
    diff_slope = 1.0 / ts_diff[0]

    # Calculate tide-corrected raw timeseries using both slopes
    corr_left = sdf_point - (tide_point / raw_slope)
    corr_right = sdf_point - (tide_point / diff_slope)
    x_corr, y_left_vals = corr_left.time.values, corr_left.values
    y_right_vals = corr_right.values

    # Standardize scatter plot aesthetics
    ts_size = 30
    sc_size = 60
    base_style = {"edgecolors": "black", "linewidths": 0.5, "alpha": 0.8, "zorder": 3}
    
    raw_style, diff_style = base_style.copy(), base_style.copy()
    raw_tide_style, diff_tide_style = base_style.copy(), base_style.copy()

    if color_by_time:
        raw_style.update({"c": sdf_point.time.dt.decimal_year.values, "cmap": "viridis"})
        diff_style.update({"c": diff_sdf.time.dt.decimal_year.values, "cmap": "viridis"})
        raw_tide_style.update(raw_style)
        diff_tide_style.update(diff_style)
    else:
        raw_style.update({"color": "steelblue"})
        diff_style.update({"color": "steelblue"})
        raw_tide_style.update({"color": "orange"})
        diff_tide_style.update({"color": "orange"})

    fig = plt.figure(figsize=(12.5, 11.5))
    gs = fig.add_gridspec(3, 6)
    
    axes = np.empty((2, 3), dtype=object)
    axes[0, 0] = fig.add_subplot(gs[0, 0:2])
    axes[0, 1] = fig.add_subplot(gs[0, 2:4])
    axes[0, 2] = fig.add_subplot(gs[0, 4:6])
    axes[1, 0] = fig.add_subplot(gs[1, 0:2])
    axes[1, 1] = fig.add_subplot(gs[1, 2:4])
    axes[1, 2] = fig.add_subplot(gs[1, 4:6])
    
    ax_corr_left = fig.add_subplot(gs[2, 0:3])
    ax_corr_right = fig.add_subplot(gs[2, 3:6])

    # Plot raw distances
    sdf_point.plot(ax=axes[0, 0], linestyle="-", color="gray", alpha=0.5)
    axes[0, 0].scatter(t_raw, y_raw, s=ts_size, **raw_style)
    axes[0, 0].set_title("Raw SDF")
    axes[0, 0].set_ylabel("Distance (m)")
    axes[0, 0].set_xlabel("time")

    # Plot raw tides
    tide_point.plot(ax=axes[0, 1], linestyle="-", color="gray", alpha=0.5)
    axes[0, 1].scatter(t_raw, x_raw, s=ts_size, **raw_tide_style)
    axes[0, 1].set_title("Raw tides")
    axes[0, 1].set_ylabel("Tide height (m)")

    # Plot raw distance vs tide scatter & trendline
    axes[0, 2].scatter(x_raw, y_raw, s=sc_size, **raw_style)
    x_bounds_raw = np.array([x_raw.min(), x_raw.max()])
    axes[0, 2].plot(
        x_bounds_raw, ts_raw[0] * x_bounds_raw + ts_raw[1], 
        color="black", linestyle="--", linewidth=2, zorder=4,
        label=f"Median slope (tan β: {raw_slope:.3f}; $R^2$: {r2_raw:.2f})"
    )
    axes[0, 2].set_xlabel("Raw tides (m)")
    axes[0, 2].set_ylabel("Raw distance (m)")
    axes[0, 2].set_title(f"Raw tide vs Raw SDF")
    axes[0, 2].legend(loc="best", framealpha=0.5)

    # Plot diff distances
    diff_sdf.plot(ax=axes[1, 0], linestyle="-", color="gray", alpha=0.5)
    axes[1, 0].scatter(t_diff, y_diff, s=ts_size, **diff_style)
    axes[1, 0].set_title("Diff SDF")
    axes[1, 0].set_ylabel("Diff distance (m)")
    axes[1, 0].set_xlabel("time")

    # Plot diff tides
    diff_tides.plot(ax=axes[1, 1], linestyle="-", color="gray", alpha=0.5)
    axes[1, 1].scatter(t_diff, x_diff, s=ts_size, **diff_tide_style)
    axes[1, 1].set_title("Diff tides")
    axes[1, 1].set_ylabel("Tide height (m)")

    # Plot diff distance vs tide scatter & trendline
    axes[1, 2].scatter(x_diff, y_diff, s=sc_size, **diff_style)
    x_bounds_diff = np.array([x_diff.min(), x_diff.max()])
    axes[1, 2].plot(
        x_bounds_diff, ts_diff[0] * x_bounds_diff + ts_diff[1], 
        color="black", linestyle="--", linewidth=2, zorder=4,
        label=f"Median slope (tan β: {diff_slope:.3f}; $R^2$: {r2_diff:.2f})"
    )
    axes[1, 2].set_xlabel("Diff tides (m)")
    axes[1, 2].set_ylabel("Diff distance (m)")
    axes[1, 2].set_title(f"Diff tide vs Diff SDF")
    axes[1, 2].legend(loc="best", framealpha=0.5)

    # Plot tide-corrected raw SDF using raw slope (left)
    corr_left.plot(ax=ax_corr_left, linestyle="-", color="gray", alpha=0.5)
    ax_corr_left.scatter(x_corr, corr_left, s=ts_size, **raw_style)
    ax_corr_left.set_title(f"Tide-corrected SDF, raw slope (σ: {corr_left.std().item():.2f} m)")
    ax_corr_left.set_ylabel("Corrected distance (m)")
    ax_corr_left.set_xlabel("time")

    # Plot tide-corrected raw SDF using diff slope (right)
    corr_right.plot(ax=ax_corr_right, linestyle="-", color="gray", alpha=0.5)
    ax_corr_right.scatter(x_corr, corr_right, s=ts_size, **raw_style)
    ax_corr_right.set_title(f"Tide-corrected SDF, diff slope (σ: {corr_right.std().item():.2f} m)")
    ax_corr_right.set_ylabel("Corrected distance (m)")
    ax_corr_right.set_xlabel("time")
    
    plt.tight_layout()
    plt.show()


def plot_pixel_regression(
    sdf_stack: xr.DataArray,
    tide_heights: xr.DataArray,
    x: float,
    y: float,
    difference: bool = False,
    mask_outliers: bool = False,
    mad_multiplier: int = 6,
    alpha: float = 0.75,
    y_limits: tuple = None,
    color_by_time: bool = False,
) -> None:
    """
    Extracts a single pixel time series, applies beach_slope_theilsen, and plots the result.

    Parameters
    ----------
    sdf_stack : xr.DataArray
        Array containing signed distance fields with dimensions of time, y, and x.
    tide_heights : xr.DataArray
        Array containing tide heights. Can be 1D (time) or 3D (time, y, x).
    x : float
        The longitude or easting of the target pixel.
    y : float
        The latitude or northing of the target pixel.
    difference : bool, optional
        If True, differences the time series before slope estimation. Default is False.
    mask_outliers : bool, optional
        If True, applies the boomerang outlier flag before differencing, 
        by default False.
    mad_multiplier : int, optional
        Multiplier for the median absolute deviation used in outlier detection,
        by default 6.
    alpha : float, optional
        Confidence degree between zero and one. Default is 0.75.
    y_limits : tuple, optional
        A tuple of (min, max) setting the y-axis boundaries for the signed distance.
        Default is None.
    color_by_time : bool, optional
        If True, colors the scatter points by decimal year to visualize temporal 
        trends within the regression. Default is False.
    """
    # Isolate single pixel
    sdf_point = sdf_stack.sel(x=x, y=y, method="nearest")
    tide_point = tide_heights.sel(x=x, y=y, method="nearest")

    # Optionally mask out outliers
    if mask_outliers:
        outlier_mask = flag_boomerang_outliers(sdf_stack=sdf_point, mad_multiplier=mad_multiplier)
        sdf_point = sdf_point.where(~outlier_mask)

    # Run Theilsen beach slope estimation, turning spatial smoothing off for one pixel
    res = beach_slope_theilsen(
        sdf_stack=sdf_point,
        tide_heights=tide_point,
        alpha=alpha,
        smooth_spatial=0,
        difference=difference,
    )

    if np.isnan(res.slope_optimal.item()):
        print("Insufficient valid data points to perform regression.")
        return

    # Extract physical beach slopes (tan beta)
    tan_beta_opt = res.slope_optimal.item()
    tan_beta_high = res.slope_high.item()
    tan_beta_low = res.slope_low.item()
    c = res.intercept.item()

    # Convert physical beach slopes to regression slopes (dx/dz) for line plotting
    alg_opt = 1.0 / tan_beta_opt
    alg_high = 1.0 / tan_beta_high
    alg_low = 1.0 / tan_beta_low

    if difference:
        sdf_point = sdf_point.diff(dim="time")
        tide_point = tide_point.diff(dim="time")

    sdf_plot, tide_plot = xr.align(
        sdf_point.dropna(dim="time", how="all"),
        tide_point.dropna(dim="time", how="all"),
        join="inner",
    )

    plt.figure(figsize=(8, 6))

    # Plot scatter points colored by time if requested
    if color_by_time:
        time_colors = sdf_plot.time.dt.decimal_year.values
        sc = plt.scatter(
            tide_plot.values,
            sdf_plot.values,
            c=time_colors,
            cmap="viridis",
            edgecolor="white",
            s=50,
            label="Observations",
        )
        cbar = plt.colorbar(sc)
        cbar.set_label("Year")
    else:
        plt.scatter(
            tide_plot.values,
            sdf_plot.values,
            c="steelblue",
            edgecolor="white",
            s=50,
            label="Observations",
        )

    # Calculate plotting bounds using the observed tidal range
    tide_min = np.nanmin(tide_plot.values)
    tide_max = np.nanmax(tide_plot.values)
    x_line = np.array([tide_min, tide_max])

    # Calculate corresponding signed distance y-values
    y_opt = alg_opt * x_line + c
    y_high = alg_high * x_line + c
    y_low = alg_low * x_line + c

    plt.fill_between(
        x_line,
        y_low,
        y_high,
        color="gray",
        alpha=0.2,
        label=f"{int(alpha*100)}% Confidence Interval",
    )
    plt.plot(
        x_line,
        y_opt,
        color="black",
        linewidth=2,
        label=f"Median Slope (tan β = {tan_beta_opt:.3f}; {tan_beta_low:.3f}-{tan_beta_high:.3f})",
    )

    # Apply custom y-limits if provided 
    if y_limits is not None:
        if isinstance(y_limits, tuple):
            plt.ylim(y_limits)
        elif y_limits == "auto":
            min_lim, max_lim = y_opt
            buffer = abs(sdf_plot - sdf_plot.median()).median() * 6
            plt.ylim(min_lim - buffer, max_lim + buffer)

    plt.title(f"Theil Sen slope estimation\nX: {x:.2f}, Y: {y:.2f}")
    plt.xlabel("Tide height (m)")
    plt.ylabel("Signed distance (m)")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.show()


def plot_pixel_dispersion(
    sdf_stack: xr.DataArray,
    tide_heights: xr.DataArray,
    x: float,
    y: float,
    difference: bool = False,
    mask_outliers: bool = False,
    mad_multiplier: int = 6,
    candidate_slopes: np.ndarray = None,
    confidence_band: float = 0.05,
    optimal_method: str = "min",
    override_slope: float | None = None,

) -> None:
    """
    Extracts a single pixel time series and visualizes the variance 
    minimization objective function, including the high-density interpolation 
    used for optimal slope extraction.
    
    Parameters
    ----------
    sdf_stack : xr.DataArray
        Array containing signed distance fields with dimensions of time, y, and x.
    tide_heights : xr.DataArray
        Array containing tide heights. Can be 1D (time) or 3D (time, y, x).
    x : float
        The longitude or easting of the target pixel.
    y : float
        The latitude or northing of the target pixel.
    difference : bool, optional
        If True, differences the time series before slope estimation. Default is True.
    mask_outliers : bool, optional
        If True, applies the boomerang outlier flag before differencing, 
        by default False.
    mad_multiplier : int, optional
        Multiplier for the median absolute deviation used in outlier detection,
        by default 6.
    candidate_slopes : np.ndarray, optional
        Array of candidate slopes to evaluate. Default is None.
    confidence_band : float, optional
        Threshold used to calculate uncertainty bounds, default 0.05.
    optimal_method : str, optional
        Method used for optimal slope extraction, default "min".
    override_slope : float, optional
        Manually override the extracted optimal slope. Default None.
    """
    # Isolate the single pixel while preserving spatial dimensions (y, x) 
    # to maintain compatibility with the pipeline slope_variance function
    sdf_point = sdf_stack.sel(x=[x], y=[y], method="nearest")

    if "x" in tide_heights.dims and "y" in tide_heights.dims:
        tide_point = tide_heights.sel(x=[x], y=[y], method="nearest")
    else:
        tide_point = tide_heights

    # Optionally mask out outliers before differencing
    if mask_outliers:
        outlier_mask = flag_boomerang_outliers(
            sdf_stack=sdf_point,
            mad_multiplier=mad_multiplier,
        )
        sdf_point = sdf_point.where(~outlier_mask)
    
    if difference:
        print("Differencing observations before slope estimation")
        sdf_point = sdf_point.diff(dim="time")
        tide_point = tide_point.diff(dim="time")

    # Align and remove missing time steps
    sdf_point, tide_point = xr.align(
        sdf_point.dropna(dim="time", how="all"),
        tide_point.dropna(dim="time", how="all"),
        join="inner",
    )
    
    if len(sdf_point.time) < 4:
        print("Insufficient valid data points to perform optimization.")
        return

    # Create a mock coastal mask covering just this single pixel
    mask = xr.DataArray(
        [[True]], 
        coords={"y": sdf_point.y, "x": sdf_point.x}, 
        dims=["y", "x"]
    )

    # Execute the production pipeline to generate the raw surface
    dispersion_surface = _slope_dispersion(
        sdf_stack=sdf_point,
        tide_heights=tide_point,
        coastal_mask=mask,
        candidate_slopes=candidate_slopes,
    )
    dispersion_surface.load()

    # Extract the final results using the production function
    slope_results = _extract_optimal_slopes(
        dispersion=dispersion_surface,
        confidence_band=confidence_band,
        optimal_method=optimal_method,
    )

    # Reuse outputs directly from slope_results
    opt_slope = slope_results.slope_optimal.item()
    slope_lower = slope_results.slope_low.item()
    slope_upper = slope_results.slope_high.item()

    # If required, overrule optimal slope
    if override_slope is not None:
        opt_slope = override_slope            

    # Isolate one-dimensional dispersion curve for plotting
    disp_1d = dispersion_surface.squeeze()

    # Calculate threshold based on the interpolated minimum to match production logic
    min_disp = disp_1d.min(dim="candidate_slope")
    threshold = min_disp * (1.0 + confidence_band)

    # Calculate final corrected timeseries for the scatter plot
    sdf_1d = sdf_point.squeeze()
    tide_1d = tide_point.squeeze()
    final_corrected_sdf = sdf_1d - (tide_1d / opt_slope)

    # Calculate Theil-Sen robust regression for both datasets
    ts_uncorr = stats.theilslopes(sdf_1d.values, tide_1d.values)
    ts_corr = stats.theilslopes(final_corrected_sdf.values, tide_1d.values)

    x_line = np.array([tide_1d.min().item(), tide_1d.max().item()])
    y_uncorr_line = ts_uncorr[0] * x_line + ts_uncorr[1]
    y_corr_line = ts_corr[0] * x_line + ts_corr[1]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot variance minimization curve
    ax1.plot(
        disp_1d.candidate_slope.values, 
        disp_1d.values, 
        color="steelblue", 
        linewidth=2.5, 
        label="Interpolated dispersions"
    )
    
    ax1.scatter(
        opt_slope, 
        min_disp, 
        color="red", 
        s=80, 
        zorder=5,
        label=f"Optimal slope (tan θ = {opt_slope:.3f})"
    )
    
    ax1.axhline(threshold, color="gray", linestyle="--", alpha=0.7, label=f"{int(confidence_band*100)}% threshold")
    ax1.axvspan(slope_lower, slope_upper, color="gray", alpha=0.15, label=f"Confidence band [{slope_lower:.3f}, {slope_upper:.3f}]")
    ax1.axvline(opt_slope, color="red", linestyle="--", alpha=0.5)
    
    ax1.set_title("Variance minimization curve")
    ax1.set_xlabel("Candidate beach slope (tan θ)")
    ax1.set_ylabel("Median Absolute Deviation (MAD)")
    ax1.legend(loc="best")

    # Plot Scatter and Correction Vectors
    ax2.scatter(
        tide_1d.values, 
        sdf_1d.values, 
        color="steelblue", 
        s=40,
        alpha=0.7,
        label="Uncorrected SDF",
        zorder=3
    )
    
    for pt_x, pt_y_orig, pt_y_corr in zip(tide_1d.values, sdf_1d.values, final_corrected_sdf.values):
        ax2.annotate(
            "",
            xy=(pt_x, pt_y_corr),
            xytext=(pt_x, pt_y_orig),
            arrowprops=dict(arrowstyle="->", color="gray", alpha=0.4, shrinkA=0, shrinkB=0),
            zorder=2
        )
        
    ax2.scatter(
        tide_1d.values, 
        final_corrected_sdf.values, 
        color="red", 
        marker="x",
        s=40,
        label="Tide-corrected SDF",
        zorder=4
    )

    # Plot the Theil-Sen regression lines
    ax2.plot(
        x_line, y_uncorr_line, 
        color="steelblue", 
        linestyle="--", 
        alpha=0.8, 
        label="Uncorrected trend",
        zorder=5
    )
    ax2.plot(
        x_line, y_corr_line, 
        color="red", 
        linestyle="--", 
        alpha=0.8, 
        label="Corrected trend",
        zorder=5
    )
    
    ax2.set_title("Tide correction")
    ax2.set_xlabel("Observed tide height (m)")
    ax2.set_ylabel("Signed distance (m)")
    ax2.legend(loc="best")

    min_lim, max_lim = y_uncorr_line 
    buffer = np.median(abs(sdf_1d.values - np.median(sdf_1d.values))) * 5
    ax2.set_ylim(min_lim - buffer, max_lim + buffer)

    fig.suptitle(f"Slope optimization, X: {x:.2f}, Y: {y:.2f}", fontsize=14)
    plt.tight_layout()
    plt.show()


def flag_boomerang_outliers(
    sdf_stack: xr.DataArray,
    mad_multiplier: float = 6.0,
    min_threshold: float = 15.0,
) -> xr.DataArray:
    """
    Identifies isolated anomalous shoreline extractions using an
    adaptive, variance-based threshold in the first-difference domain.

    Parameters
    ----------
    sdf_stack : xr.DataArray
        Array of Signed Distance Fields with a time dimension.
    mad_multiplier : float, optional
        Scalar applied to the Median Absolute Deviation to define
        the outlier cutoff.
    min_threshold : float, optional
        Absolute minimum threshold in meters to prevent over-flagging
        in highly stable, micro-tidal environments.

    Returns
    -------
    xr.DataArray
        Boolean mask where True indicates a detected outlier.
    """
    diff_forward = sdf_stack.diff(dim="time")

    # Shift array backward to align consecutive differences
    diff_subsequent = diff_forward.shift(time=-1, fill_value=0)

    # Calculate robust pixel-wise dispersion of shoreline differences
    median_diff = diff_forward.median(dim="time", skipna=True)
    mad_diff = np.abs(diff_forward - median_diff).median(dim="time", skipna=True)

    # Establish adaptive threshold with a hard floor
    adaptive_threshold = mad_diff * mad_multiplier
    dynamic_limit = adaptive_threshold.clip(min=min_threshold)

    # Enforce that jumps exceed the local dynamic limit
    is_massive_forward = np.abs(diff_forward) > dynamic_limit
    is_massive_subsequent = np.abs(diff_subsequent) > dynamic_limit

    # Confirm the shoreline snapped back
    reverses_direction = np.sign(diff_forward) != np.sign(diff_subsequent)

    # Flag as outlier
    is_outlier = is_massive_forward & is_massive_subsequent & reverses_direction

    # Reindex mask back to original coordinate alignment
    return is_outlier.reindex_like(sdf_stack, fill_value=False)