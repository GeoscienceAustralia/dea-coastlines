import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

from coastlines.timeseries import _block_theil_sen, slope_variance, extract_optimal_slopes


def plot_pixel_regression(
    sdf_stack: xr.DataArray,
    tide_heights: xr.DataArray,
    x_coord: float,
    y_coord: float,
    alpha: float = 0.75,
    y_limits: tuple = None,
    color_by_time: bool = False,
) -> None:
    """
    Extracts a single pixel time series, applies the block-level
    Theil-Sen regression used in the broader pipeline, and plots the result.

    Parameters
    ----------
    sdf_stack : xr.DataArray
        Array containing signed distance fields with dimensions of time, y, and x.
    tide_heights : xr.DataArray
        Array containing tide heights. Can be 1D (time) or 3D (time, y, x).
    x_coord : float
        The longitude or easting of the pixel to plot.
    y_coord : float
        The latitude or northing of the pixel to plot.
    alpha : float, optional
        Confidence degree between zero and one. Default is 0.75.
    y_limits : tuple, optional
        A tuple of (min, max) setting the y-axis boundaries for the signed distance.
        Default is None.
    color_by_time : bool, optional
        If True, colors the scatter points by decimal year to visualize temporal 
        trends within the regression. Default is False.
    """
    # Isolate the spatial point for the SDF using nearest neighbour interpolation
    sdf_point = sdf_stack.sel(x=x_coord, y=y_coord, method="nearest")

    # Safely isolate the tide point if the tide array is spatially varying
    if "x" in tide_heights.dims and "y" in tide_heights.dims:
        tide_point_raw = tide_heights.sel(x=x_coord, y=y_coord, method="nearest")
    else:
        tide_point_raw = tide_heights

    # Align and drop obscured timesteps to precisely mirror the wrapper function
    sdf_point, tide_point = xr.align(
        sdf_point.dropna(dim="time", how="all"),
        tide_point_raw.dropna(dim="time", how="all"),
        join="inner",
    )

    # Ensure there is sufficient data to run a meaningful regression
    valid_points = ~np.isnan(sdf_point.values) & ~np.isnan(tide_point.values)
    if np.sum(valid_points) < 4:
        print("Insufficient valid data points to perform regression.")
        return

    # Restructure the flat arrays into pseudo-spatial blocks to feed
    # directly into the existing block helper function
    y_chunk = sdf_point.values[np.newaxis, np.newaxis, :]
    x_chunk = tide_point.values

    # Execute the regression logic used in the parallel pipeline
    opt, high, low, intercept = _block_theil_sen(y_chunk, x_chunk, alpha=alpha)

    # Unpack the scalar results 
    alg_opt = opt[0, 0]
    alg_high = high[0, 0]
    alg_low = low[0, 0]
    c = intercept[0, 0]

    # Initialize the plot layout
    plt.figure(figsize=(8, 6))
    
    # Plot scatter points colored by time if requested
    if color_by_time:
        # Extract decimal years for the colormap
        time_colors = sdf_point.time.dt.decimal_year.values
        sc = plt.scatter(
            tide_point.values,
            sdf_point.values,
            c=time_colors,
            cmap="viridis",
            edgecolor="white",
            s=50,
            label="Observations",
        )
        # Add a colorbar to indicate the time gradient
        cbar = plt.colorbar(sc)
        cbar.set_label("Year")
    else:
        # Fall back to a solid color if time coloring is disabled
        plt.scatter(
            tide_point.values,
            sdf_point.values,
            c="steelblue",
            edgecolor="white",
            s=50,
            label="Observations",
        )
  
    # Calculate plotting bounds using the observed tidal range
    tide_min = np.nanmin(tide_point.values)
    tide_max = np.nanmax(tide_point.values)
    x_line = np.array([tide_min, tide_max])

    # Calculate corresponding signed distance y-values
    y_opt = alg_opt * x_line + c
    y_high = alg_high * x_line + c
    y_low = alg_low * x_line + c

    # Get physical slope
    slope = 1.0 / alg_opt
    slope_min = 1.0 / alg_low
    slope_max = 1.0 / alg_high

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
        label=f"Median Slope (tan β = {slope:.2f}; {slope_min:.2f}-{slope_max:.2f})",
    )

    # Apply custom y-limits if provided 
    if y_limits is not None:
        if isinstance(y_limits, tuple):
            plt.ylim(y_limits)
        if y_limits == "auto":
            min_lim, max_lim = y_opt
            buffer = abs(sdf_point - sdf_point.median()).median() * 5
            plt.ylim(min_lim - buffer, max_lim + buffer)

    # Finalize plot aesthetics
    plt.title(f"Theil Sen slope estimation\nX: {x_coord:.2f}, Y: {y_coord:.2f}")
    plt.xlabel("Tide height (m)")
    plt.ylabel("Signed distance (m)")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.show()


def plot_slope_variance(
    sdf_stack: xr.DataArray,
    tide_heights: xr.DataArray,
    x_coord: float,
    y_coord: float,
    candidate_slopes: np.ndarray = None,
    confidence_band: float = 0.05,
    interp_density: int = 0,
    interp_method: str = "linear",
) -> None:
    """
    Extracts a single pixel time series and visualizes the variance 
    minimization objective function, including the high-density interpolation 
    used for optimal slope extraction.
    """
    # Isolate the single pixel while preserving spatial dimensions (y, x) 
    sdf_point = sdf_stack.sel(x=[x_coord], y=[y_coord], method="nearest")

    if "x" in tide_heights.dims and "y" in tide_heights.dims:
        tide_point = tide_heights.sel(x=[x_coord], y=[y_coord], method="nearest")
    else:
        tide_point = tide_heights

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
    dispersion_surface = slope_variance(
        sdf_stack=sdf_point,
        tide_heights=tide_point,
        coastal_mask=mask,
        candidate_slopes=candidate_slopes
    )
    dispersion_surface.load()

    # Extract the final results using the production function
    slope_results = extract_optimal_slopes(
        dispersion=dispersion_surface,
        confidence_band=confidence_band,
        smooth_window=0,  # Disable spatial smoothing for a single point
        interp_density=interp_density,
        interp_method=interp_method,
    )

    # Reuse outputs directly from slope_results
    opt_slope = slope_results.slope_optimal.item()
    slope_lower = slope_results.slope_low.item()
    slope_upper = slope_results.slope_high.item()

    # Isolate the one-dimensional raw dispersion curve for plotting
    raw_disp = dispersion_surface.squeeze()

    # Recreate the high-density interpolation to match the extraction process exactly
    if interp_density > 0:
        min_slope_val = float(raw_disp.candidate_slope.min())
        max_slope_val = float(raw_disp.candidate_slope.max())
        dense_slopes = np.linspace(min_slope_val, max_slope_val, interp_density, dtype=np.float32)
        interp_disp = raw_disp.interp(candidate_slope=dense_slopes, method=interp_method)
    else: 
        interp_disp = raw_disp

    # Calculate threshold based on the interpolated minimum to match production logic
    min_disp = interp_disp.min().item()
    # threshold = min_disp * (1.0 + confidence_band)
    # min_disp = interp_disp.min().item()
    # max_disp = interp_disp.max().item()
    # threshold = min_disp + (confidence_band * (max_disp - min_disp))
    threshold = interp_disp.quantile(q=confidence_band)

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

    # Set up the visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Panel A: Objective Function Curve
    ax1.plot(
        raw_disp.candidate_slope.values, 
        raw_disp.values,
        color="gray", 
        marker="o",
        linewidth=0,
        alpha=0.5,
        label="Raw dispersions"
    )
    
    ax1.plot(
        interp_disp.candidate_slope.values, 
        interp_disp.values, 
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

    # Panel B: Scatter & Correction Vectors
    ax2.scatter(
        tide_1d.values, 
        sdf_1d.values, 
        color="steelblue", 
        s=40,
        alpha=0.7,
        label="Uncorrected SDF",
        zorder=3
    )
    
    for x, y_orig, y_corr in zip(tide_1d.values, sdf_1d.values, final_corrected_sdf.values):
        ax2.annotate(
            "",
            xy=(x, y_corr),
            xytext=(x, y_orig),
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
        label=f"Uncorrected trend",
        zorder=5
    )
    ax2.plot(
        x_line, y_corr_line, 
        color="red", 
        linestyle="--", 
        alpha=0.8, 
        label=f"Corrected trend",
        zorder=5
    )
    
    ax2.set_title("Tide correction")
    ax2.set_xlabel("Observed tide height (m)")
    ax2.set_ylabel("Signed distance (m)")
    ax2.legend(loc="best")

    min_lim, max_lim = y_uncorr_line 
    buffer = np.median(abs(sdf_1d.values - np.median(sdf_1d.values))) * 5
    ax2.set_ylim(min_lim - buffer, max_lim + buffer)

    fig.suptitle(f"Slope optimization, X: {x_coord:.2f}, Y: {y_coord:.2f}", fontsize=14)
    plt.tight_layout()
    plt.show()
