import numpy as np
import xarray as xr
from scipy.spatial import cKDTree
from skimage.measure import find_contours
from scipy.stats import theilslopes
import matplotlib.pyplot as plt
from scipy import stats


def exact_subpixel_sdf(mndwi_da, threshold=0.0, resolution=10.0, min_vertices=10):
    """
    Calculates an exact sub-pixel Signed Distance Field by extracting
    marching squares contours and querying a KDTree. Optimised for speed
    using valid-pixel masking and multi-threading.

    Parameters
    ----------
    mndwi_da : xr.DataArray
        A 2D DataArray containing the calculated MNDWI values.
    threshold : float, optional
        The zero-crossing threshold. Default is 0.0.
    resolution : float, optional
        The pixel size in metres. Default is 10.0.
    min_vertices : int, optional
        The minimum number of vertices required to retain a contour string,
        used to filter out small noisey loops. Default is 10.

    Returns
    -------
    xr.DataArray
        A 2D DataArray containing the sub-pixel SDF in metres, or NaNs if
        no shoreline interface exists.
    """
    # Extract values
    mndwi_da = mndwi_da.squeeze()
    vals = mndwi_da.values

    # Extract subpixel precision contour vertices using marching squares
    contours = find_contours(vals, level=threshold)
    contours = [c for c in contours if len(c) > min_vertices]

    # Handle scenes with no water/land interface
    if not contours:
        nan_array = xr.full_like(mndwi_da, np.nan, dtype=np.float32)
        nan_array.name = "sdf_exact"
        return nan_array

    # Safely stack contours now that we know the list is not empty
    contour_pts = np.vstack(contours)

    # Filter out any non-finite values (NaN or inf)
    valid_pts_mask = np.isfinite(contour_pts).all(axis=1)
    contour_pts = contour_pts[valid_pts_mask]

    # Build a KDTree from the shoreline points
    tree = cKDTree(contour_pts)

    # Create a boolean mask of valid (non-NaN) pixels, and
    # extract integer grid coordinates for valid pixels only
    valid_mask = ~np.isnan(vals)
    valid_y, valid_x = np.nonzero(valid_mask)
    valid_grid_pts = np.c_[valid_y, valid_x]

    # Query the tree only for valid coordinates using multithreading
    distances, _ = tree.query(valid_grid_pts, workers=-1)

    # Initialize the empty array as float32 to save memory
    sdf_vals = np.full(vals.shape, np.nan, dtype=np.float32)

    # Convert to physical distances
    abs_dist_metres = distances * resolution

    # Convert to signed distances based on original MNDWI signs
    # Positive for land (MNDWI >= threshold), negative for water (MNDWI < threshold)
    valid_signs = np.where(vals[valid_mask] >= threshold, 1.0, -1.0)
    sdf_vals[valid_mask] = abs_dist_metres * valid_signs

    return xr.DataArray(
        sdf_vals, coords=mndwi_da.coords, dims=mndwi_da.dims, name="sdf_exact"
    )

import numpy as np
import xarray as xr
import skfmm

def subpixel_sdf_fmm(da, threshold=0.0, resolution=10.0):
    """
    Calculates an exact sub-pixel Signed Distance Field (SDF)
    using the Fast Marching Method (scikit-fmm).

    Parameters
    ----------
    da : xr.DataArray
        A 2D DataArray to convert to a sub-pixel SDF.
    threshold : float, optional
        The zero-crossing threshold. Default is 0.0.
    resolution : float, optional
        The pixel size in metres. Default is 10.0.

    Returns
    -------
    xr.DataArray
        A 2D DataArray containing the sub-pixel SDF in metres,
        or NaNs if zero-crossing interface exists.
    """
    # Extract as numpy arrays
    da = da.squeeze()
    vals = da.values

    # Subtract thresholds to support 0-threshold extraction
    phi = vals - threshold

    # Create valid data mask to propogate FMM around NaNs
    valid_mask = ~np.isnan(phi)

    # Return empty array if no 0 interface exists in data
    if not (np.any(phi[valid_mask] > 0) and np.any(phi[valid_mask] <= 0)):
        nan_array = xr.full_like(da, np.nan, dtype=np.float32)
        nan_array.name = "sdf_exact"
        return nan_array

    # Convert to a masked array for NaN handing
    phi_masked = np.ma.masked_where(~valid_mask, phi)

    # Calculate signed distance field in physical metre units
    sdf_masked = skfmm.distance(phi_masked, dx=resolution)

    # Re-apply NaNs and cast to float32 to reduce memory
    sdf_vals = sdf_masked.filled(np.nan).astype(np.float32)

    return xr.DataArray(
        sdf_vals, coords=da.coords, dims=da.dims, name="sdf_exact"
    )


def _smooth_spatial_array(
    da: xr.DataArray | xr.Dataset, method: str = "median", window_size: int = 3
) -> xr.DataArray | xr.Dataset:
    """
    Applies a spatial median filter to 2D arrys to reduce noise
    while preserving sharp morphological edges.

    Parameters
    ----------
    da : xr.DataArray | xr.Dataset
        An xr.DataArray containing 2D spatial data to smooth, or a
        xr.Dataset containing multiple 2D variables to smooth.
    method : str, optional
        Whether to use "median" or "mean" smoothing.
    window_size : int, optional
        The size of the moving window. Must be an odd integer. Default is 3.

    Returns
    -------
    xr.DataArray | xr.Dataset
        Smoothed output arrays.
    """
    if window_size % 2 == 0:
        raise ValueError("Window size must be an odd integer.")

    if method == "median":
        smoothed = (
            da.rolling(y=window_size, x=window_size, center=True)
            .construct(y="window_y", x="window_x")
            .median(dim=["window_y", "window_x"], skipna=True)
        )
    elif method == "mean":
        smoothed = (
            da.rolling(y=window_size, x=window_size, center=True)
            .construct(y="window_y", x="window_x")
            .mean(dim=["window_y", "window_x"], skipna=True)
        )

    return smoothed


def _block_theil_sen(y_chunk: np.ndarray, x_chunk: np.ndarray, alpha: float) -> tuple:
    """
    Block-level helper function to process generic Theil-Sen regression over
    spatial dimensions, returning raw gradients and intercepts.
    """
    # xarray apply_ufunc moves the core dimension ('time') to the last axis
    # If the input was (time, y, x), it enters this function as (y, x, time)
    spatial_shape = y_chunk.shape[:-1]
    time_dim = y_chunk.shape[-1]

    # Flatten all spatial dimensions to a single axis for easy iteration
    y_flat = y_chunk.reshape(-1, time_dim)

    # Handle either 1D (time,) or (N_pixels, time) inputs
    if x_chunk.ndim > 1:
        x_flat = x_chunk.reshape(-1, time_dim)
    else:
        x_flat = x_chunk

    # Initialize flat output arrays using float32 to save memory
    n_pixels = y_flat.shape[0]
    m_opt = np.full(n_pixels, np.nan, dtype=np.float32)
    m_l = np.full(n_pixels, np.nan, dtype=np.float32)
    m_h = np.full(n_pixels, np.nan, dtype=np.float32)
    intercept = np.full(n_pixels, np.nan, dtype=np.float32)

    # Create a mask of valid pixels requiring at least 4 valid time steps
    valid_counts = np.sum(~np.isnan(y_flat), axis=1)
    valid_pixels = np.nonzero(valid_counts >= 4)[0]

    # Halt execution instantly if the block contains no valid data
    if valid_pixels.size == 0:
        return (
            m_opt.reshape(spatial_shape),
            m_l.reshape(spatial_shape),
            m_h.reshape(spatial_shape),
            intercept.reshape(spatial_shape),
        )

    # Iterate strictly over valid pixels
    for idx in valid_pixels:
        y_ts = y_flat[idx]
        x_ts = x_flat[idx] if x_chunk.ndim > 1 else x_flat

        valid = ~np.isnan(x_ts) & ~np.isnan(y_ts)

        # Retrieve raw gradient and intercept
        # Apply alpha to set confidence bounds (e.g. 0.90 == 90% confidence)
        m, c, m_low, m_high = theilslopes(y_ts[valid], x_ts[valid], alpha=alpha)

        m_opt[idx] = m
        m_l[idx] = m_low
        m_h[idx] = m_high
        intercept[idx] = c

    # Reshape the 1D results back to the original 2D spatial dimensions
    return (
        m_opt.reshape(spatial_shape),
        m_l.reshape(spatial_shape),
        m_h.reshape(spatial_shape),
        intercept.reshape(spatial_shape),
    )


def beach_slope_theilsen(
    sdf_stack: xr.DataArray,
    tide_heights: xr.DataArray,
    alpha: float = 0.75,
    smooth_spatial: int = 3,
    smooth_method: str = "mean",
    difference: bool = False,
) -> xr.Dataset:
    """
    Calculates the robust beach slope, intercept, and confidence bounds.

    Parameters
    ----------
    sdf_stack : xr.DataArray
        A 3D array (time, y, x) of signed distance fields.
    tide_heights : xr.DataArray
        A 1D (time) or 3D  (time, y, x) array of tide heights.
    alpha : float, optional
        Confidence degree between 0 and 1. Default is 0.75 for a 75%
        confidence interval.
    smooth_spatial : int, optional
        The size of the spatial filter to apply to the outputs to reduce
        spatial noise. Set to zero to disable smoothing. Default is 3.
    smooth_method : str, optional
        Method to use for spatial smoothing; supports "mean" and "median".
        Default is "mean".
    difference : bool, optional
        Whether to apply differencing to SDF and tide data along the
        "time" axis prior to slope generation. This serves as a
        high-pass filter, retaining high frequency tidal change,
        but supressing lower frequency long-term coastal change.

    Returns
    -------
    xr.Dataset
        A dataset containing 2D maps: median slope, steepest probable slope,
        flattest probable slope, regression intercept, and slope uncertainty.
    """
    # Drop any redundant all-NaN timesteps to improve run-time
    sdf_stack, tide_heights = xr.align(
        sdf_stack.dropna(dim="time", how="all"),
        tide_heights.dropna(dim="time", how="all"),
        join="inner",
    )

    # Optionally calculate differences between consecutive observations,
    # rather than operate on absolute shoreline positions. This
    # supresses low frequency long-term coastal change, but retains
    # higher frequency tidal change.
    if difference:
        print("Differencing observations before slope estimation")
        sdf_stack = sdf_stack.diff(dim="time")
        tide_heights = tide_heights.diff(dim="time")

    # Run Theil Sen gradient estimation in parallel
    raw_m, raw_m_l, raw_m_h, intercept = xr.apply_ufunc(
        _block_theil_sen,
        sdf_stack,
        tide_heights,
        kwargs={"alpha": alpha},
        input_core_dims=[["time"], ["time"]],
        output_core_dims=[[], [], [], []],
        vectorize=False,
        dask="parallelized",
        output_dtypes=[np.float32, np.float32, np.float32, np.float32],
    )

    # Invert raw gradients to physical beach slopes.
    # The lowest gradient (m_l) corresponds to the steepest slope (slope_high).
    slope_optimal = xr.where(raw_m != 0, 1.0 / raw_m, np.nan)
    slope_high = xr.where(raw_m_l != 0, 1.0 / raw_m_l, np.nan)
    slope_low = xr.where(raw_m_h != 0, 1.0 / raw_m_h, np.nan)

    # Package raw geometric arrays into an xarray Dataset
    ds_out = xr.Dataset(
        {
            "slope_optimal": slope_optimal,
            "slope_high": slope_high,
            "slope_low": slope_low,
            "intercept": intercept,
        }
    )

    # Apply spatial filtering to all dataset variables if requested
    if smooth_spatial > 0:
        ds_out = _smooth_spatial_array(
            ds_out,
            method=smooth_method,
            window_size=smooth_spatial,
        )

    # Calculate final uncertainty safely using the smoothed upper and lower bounds
    ds_out["slope_uncertainty"] = np.abs(ds_out.slope_low - ds_out.slope_high)

    return ds_out


def _slope_dispersion(
    sdf_stack: xr.DataArray,
    tide_heights: xr.DataArray,
    coastal_mask: xr.DataArray,
    candidate_slopes: np.ndarray = None,
) -> xr.DataArray:
    """
    Calculates shoreline position dispersion across a range of candidate
    beach slopes. Dispersion is measured by Median Absolute Deviation (MAD).

    Parameters
    ----------
    sdf_stack : xr.DataArray
        Array of uncorrected Signed Distance Fields with dimensions time, y, and x.
    tide_heights : xr.DataArray
        Array of tide heights with dimensions time, y, and x.
    coastal_mask : xr.DataArray
        Boolean mask array to restrict processing to valid coastal areas.
    candidate_slopes : np.ndarray, optional
        Array of physical slopes to test. If none are provided, a non-linear
        distribution focusing on flatter slopes is used.

    Returns
    -------
    xr.DataArray
        A multidimensional array containing the computed MAD dispersion
        metric for each candidate slope, masked to the coastal zone.
    """
    if candidate_slopes is None:
        candidate_slopes = np.geomspace(0.005, 0.3, num=40)

    test_slopes = xr.DataArray(
        candidate_slopes,
        dims=["candidate_slope"],
        coords={"candidate_slope": candidate_slopes},
    )

    # Broadcast shifts across all candidate slopes simultaneously
    horizontal_shifts = tide_heights / test_slopes
    corrected_sdfs = sdf_stack - horizontal_shifts

    # Calculate robust Median Absolute Deviation directly on the multidimensional stack
    median_val = corrected_sdfs.median(dim="time")
    dispersion = np.abs(corrected_sdfs - median_val).median(dim="time")

    # Restrict output strictly to valid data areas to save memory downstream
    return dispersion.where(coastal_mask)


def _extract_optimal_slopes(
    dispersion: xr.DataArray,
    confidence_band: float = 0.05,
    optimal_method: str = "min",
    smooth_rolling: int = 3,
) -> xr.Dataset:
    """
    Extracts the optimal slope and calculates confidence intervals from a
    computed dispersion surface, applying spatial smoothing to the outputs.

    Parameters
    ----------
    dispersion : xr.DataArray
        The pre-computed dispersion surface. This array should ideally be loaded
        into memory prior to calling this function.
    confidence_band : float, optional
        Threshold used to calculate uncertainty bounds, represented as a
        percentage above the minimum (i.e. optimal) dispersion.
    optimal_method : str, optional
        Method used to identify the optimal slope value that minimises
        dispersion after correcting for slope. Supports "min" (the input
        candidate slope that minimises variance), "mean" (mean of all slopes
        within the confidence bounds), "gmean" (geometric mean of all slopes
        within confidence bounds), and "subgrid" (a more precise version of
        "min" that interpolates an optimal slope at better than the resolution
        of the input candidate slopes).
    smooth_rolling : int, optional
        Size of the rolling mean used to smooth dispersion values along the
        "candidate_slope" axis prior to optimal  slope extraction. Set to 0
        to apply no smoothing. Default is 9.

    Returns
    -------
    xr.Dataset
        Dataset containing the optimal slope, lower bound, upper bound,
        and the absolute slope uncertainty.
    """
    # Identify valid pixels to restore NaNs later
    valid_pixels = dispersion.notnull().any(dim="candidate_slope")

    # Compute rolling mean
    if smooth_rolling > 0:
        print("Applying rolling mean")
        dispersion = dispersion.rolling(
            candidate_slope=smooth_rolling, min_periods=1, center=True
        ).mean()

    # Identify minimum dispersion and condifence bounds around it
    min_dispersion = dispersion.min(dim="candidate_slope")
    threshold = min_dispersion * (1.0 + confidence_band)
    valid_candidate_slopes = dispersion.candidate_slope.where(dispersion <= threshold)
    slope_lower = valid_candidate_slopes.min(dim="candidate_slope")
    slope_upper = valid_candidate_slopes.max(dim="candidate_slope")

    # Fill NaNs with infinity to prevent idxmin/argmin from crashing on all-NaN slices
    safe_dispersion = dispersion.fillna(np.inf)

    # Extract the optimal slope that minimises absolute dispersion
    if optimal_method == "min":
        optimal_slopes = safe_dispersion.idxmin(dim="candidate_slope")
    elif optimal_method == "mean":
        optimal_slopes = valid_candidate_slopes.mean(dim="candidate_slope")
    elif optimal_method == "gmean":
        optimal_slopes = np.exp(
            np.log(valid_candidate_slopes).mean(dim="candidate_slope")
        )
    elif optimal_method == "subgrid":
        # Get the discrete slopes as our safe fallback
        discrete_slopes = safe_dispersion.idxmin(dim="candidate_slope")
        min_idx = safe_dispersion.argmin(dim="candidate_slope")

        # Clamp indices to prevent .isel() from crashing (we mask invalid pixels later)
        max_idx = len(dispersion.candidate_slope) - 1
        safe_idx = min_idx.clip(1, max_idx - 1)

        # Extract coordinates
        m_mid = dispersion.candidate_slope.isel(candidate_slope=safe_idx)
        m_low = dispersion.candidate_slope.isel(candidate_slope=safe_idx - 1)
        m_high = dispersion.candidate_slope.isel(candidate_slope=safe_idx + 1)

        # Extract dispersion values (using original array so true NaNs are preserved)
        y_mid = dispersion.isel(candidate_slope=safe_idx)
        y_low = dispersion.isel(candidate_slope=safe_idx - 1)
        y_high = dispersion.isel(candidate_slope=safe_idx + 1)

        # Identify valid observations that are not on boundary and with three valid neighbors
        valid_subgrid = (
            (min_idx > 0) & (min_idx < max_idx) & y_low.notnull() & y_high.notnull()
        )

        # Calculate analytical vertex offset
        delta_low = m_mid - m_low
        delta_high = m_mid - m_high
        num = (delta_low**2) * (y_mid - y_high) - (delta_high**2) * (y_mid - y_low)
        den = delta_low * (y_mid - y_high) - delta_high * (y_mid - y_low)
        safe_den = den.where(den != 0, np.nan)
        vertex_offset = 0.5 * (num / safe_den).fillna(0.0)

        # Bound to physical neighborhood to prevent math explosion on noise spikes
        vertex_offset = vertex_offset.clip(min=delta_high, max=delta_low)
        subgrid_slopes = m_mid - vertex_offset

        # Merge: Use subgrid if valid, otherwise fallback to the exact discrete slope
        optimal_slopes = xr.where(valid_subgrid, subgrid_slopes, discrete_slopes)

        # Final bounds clamp
        min_bound = float(dispersion.candidate_slope.min())
        max_bound = float(dispersion.candidate_slope.max())
        optimal_slopes = optimal_slopes.clip(min_bound, max_bound)

    # Restore NaNs to pixels that contained only NoData
    optimal_slopes = optimal_slopes.where(valid_pixels)

    # Package findings into a structured dataset
    var_prob = xr.Dataset(
        {
            "slope_optimal": optimal_slopes,
            "slope_low": slope_lower,
            "slope_high": slope_upper,
        }
    )

    return var_prob


def beach_slope_dispersion(
    sdf_stack: xr.DataArray,
    tide_heights: xr.DataArray,
    coastal_mask: xr.DataArray,
    candidate_slopes: np.ndarray = None,
    confidence_band: float = 0.05,
    optimal_method: str = "subgrid",
    smooth_rolling: int = 3,
    smooth_spatial: int = 9,
    smooth_method: str = "mean",
    difference: bool = False,
) -> xr.Dataset:
    """
    Calculates the robust beach slope and confidence bounds using a 
    variance dispersion grid-search approach.

    This function wraps `slope_dispersion` and `extract_optimal_slopes` 
    to provide an API consistent with `beach_slope_theilsen`.

    Parameters
    ----------
    sdf_stack : xr.DataArray
        A 3D array (time, y, x) of signed distance fields.
    tide_heights : xr.DataArray
        A 1D (time) or 3D (time, y, x) array of tide heights.
    coastal_mask : xr.DataArray
        Boolean mask array to restrict processing to valid coastal areas.
    candidate_slopes : np.ndarray, optional
        Array of physical slopes to test. If none are provided, a non-linear
        distribution focusing on flatter slopes is used.
    confidence_band : float, optional
        Threshold used to calculate uncertainty bounds, represented as a
        percentage above the minimum (i.e. optimal) dispersion. Default is 0.05.
    optimal_method : str, optional
        Method used to identify the optimal slope value. Default is "subgrid".
    smooth_rolling : int, optional
        Size of the rolling mean used to smooth dispersion values along the
        "candidate_slope" axis prior to optimal  slope extraction. Set to 0
        to apply no smoothing. Default is 9.
    smooth_spatial : int, optional
        The size of the spatial filter to apply to the outputs to reduce
        spatial noise. Set to zero to disable smoothing. Default is 3.
    smooth_method : str, optional
        Method to use for spatial smoothing; supports "mean" and "median".
        Default is "mean".
    difference : bool, optional
        Whether to apply differencing to SDF and tide data along the
        "time" axis prior to slope generation. Default is False.

    Returns
    -------
    xr.Dataset
        A dataset containing the optimal slope, lower bound, upper bound,
        and the absolute slope uncertainty.
    """
    # Drop any redundant all-NaN timesteps to improve run-time
    sdf_stack, tide_heights = xr.align(
        sdf_stack.dropna(dim="time", how="all"),
        tide_heights.dropna(dim="time", how="all"),
        join="inner",
    )

    # Optionally calculate differences between consecutive observations,
    # rather than operate on absolute shoreline positions. This
    # supresses low frequency long-term coastal change, but retains
    # higher frequency tidal change.
    if difference:
        print("Differencing observations before slope estimation")
        sdf_stack = sdf_stack.diff(dim="time")
        tide_heights = tide_heights.diff(dim="time")

    # Generate multidimensional dispersion array
    lazy_dispersion = _slope_dispersion(
        sdf_stack=sdf_stack,
        tide_heights=tide_heights,
        coastal_mask=coastal_mask,
        candidate_slopes=candidate_slopes,
    )

    # Load dispersion grid into memory
    computed_dispersion = lazy_dispersion.compute()

    # Extract optimal slopes and confidence bounds
    ds_out = _extract_optimal_slopes(
        dispersion=computed_dispersion,
        confidence_band=confidence_band,
        optimal_method=optimal_method,
        smooth_rolling=smooth_rolling,
    )
    
    # Apply spatial filtering to all dataset variables if requested
    if smooth_spatial > 0:
        ds_out = _smooth_spatial_array(
            ds_out,
            method=smooth_method,
            window_size=smooth_spatial,
        )

    # Calculate final uncertainty safely using the smoothed bounds
    ds_out["slope_uncertainty"] = np.abs(ds_out.slope_low - ds_out.slope_high)

    return ds_out

    

def tide_correction(
    sdf_stack: xr.DataArray,
    tide_heights: xr.DataArray,
    slope_map: xr.DataArray,
    tide_datum: float = 0.0,
) -> xr.DataArray:
    """
    Applies a physical beach slope to a temporal stack of Signed Distance Fields
    to geometrically shift instantaneous shorelines to a specified tide datum.

    Parameters
    ----------
    sdf_stack : xr.DataArray
        A 3D array (time, y, x) of uncorrected SDFs.
    tide_heights : xr.DataArray
        A 1D (time) or 3D (time, y, x) array  of tide heights in meters.
    slope_map : xr.DataArray
        A 2D array (y, x) of physical beach slopes (tan θ).
    tide_datum : float, optional
        The vertical tide datum in meters to shift the shorelines toward.
        Default is 0.0.

    Returns
    -------
    xr.DataArray
        A 3D array (time, y, x) of tide-corrected SDFs.
    """
    # Calculate vertical difference between observed tide & target tide datum
    tide_difference = tide_heights - tide_datum

    # Calculate the horizontal shift for every pixel at every timestep
    # Resulting dims: (time, y, x) due to xarray broadcasting
    horizontal_offsets = tide_difference / slope_map

    # Subtract the offset from the raw distance fields
    corrected_sdf_stack = sdf_stack - horizontal_offsets

    # If pixels are NaN (e.g. no valid slope could be calculated),
    # fill them with the original uncorrected SDF stack
    # TODO: keep record of this and flag in output shoreline metadata
    corrected_sdf_stack = corrected_sdf_stack.fillna(sdf_stack)

    return corrected_sdf_stack
