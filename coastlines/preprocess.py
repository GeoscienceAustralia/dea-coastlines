import xarray as xr
import numpy as np
from skimage.measure import label, regionprops
from rasterio.features import sieve
from odc.algo import mask_cleanup


def _process_inland_mask_2d(
    water_slice: np.ndarray,
    ocean_slice: np.ndarray,
    connectivity: int,
) -> np.ndarray:
    """
    Helper function applying labelling and region properties to a single
    2D spatial array to isolate inland waterbodies.
    """
    # Break spatial array into unique, discrete regions/blobs
    blobs = label(water_slice, connectivity=connectivity)

    # Extract labels that do not overlap with the static ocean seed
    valid_labels = [
        i.label
        for i in regionprops(blobs, intensity_image=ocean_slice)
        if i.max_intensity == 0
    ]

    # Create the mask from valid inland regions
    inland_mask = np.isin(blobs, valid_labels)

    return inland_mask


def inland_water_masking(
    water_mask: xr.DataArray,
    ocean_mask: xr.DataArray,
    connectivity: int = 1,
    dilation: int = 0,
) -> xr.DataArray:
    """
    Identifies inland water by selecting regions of water that do not overlap
    with ocean pixels. This region can be optionally dilated.

    Parameters
    ----------
    water_mask : xarray.DataArray
        An array containing boolean values where 1 == water and 0 == land.
        Can be 2D or 3D (time, y, x).
    ocean_mask : xarray.DataArray
        A supplementary static boolean dataset used to separate ocean waters
        from other inland water. The array should contain values of 1
        for high certainty ocean pixels, and 0 for all other pixels.
    connectivity : int, optional
        Passed to the 'connectivity' parameter of `skimage.measure.label`. Default is 1.
    dilation : int, optional
        Number of iterations to dilate the final inland water mask. Default is 0.

    Returns
    -------
    xarray.DataArray
        An array containing the mask consisting of identified inland water
        pixels as True.
    """
    # Map the 2D helper across the spatial dimensions of the dataset
    inland_masked = xr.apply_ufunc(
        _process_inland_mask_2d,
        water_mask,
        ocean_mask,
        kwargs={"connectivity": connectivity},
        input_core_dims=[["y", "x"], ["y", "x"]],
        output_core_dims=[["y", "x"]],
        vectorize=True,
        dask="parallelized",
        output_dtypes=[bool],
    )

    inland_masked.name = "inland_water_mask"

    # Optionally dilate the resulting mask
    if dilation > 0:
        inland_masked = mask_cleanup(
            inland_masked, mask_filters=[("dilation", dilation)]
        )

    return inland_masked.astype(bool)


def create_coastal_buffer(mndwi, inland_water_mask, dilation_pixels=10, sieve_size=9):
    """Calculate a robust coastal buffer spanning the intertidal zone.

    Parameters
    ----------
    mndwi : xarray.DataArray
        Time series of MNDWI images.
    inland_water_mask : xarray.DataArray
        Boolean mask where True indicates inland water bodies.
    dilation_pixels : int, default 10
        Number of pixels to dilate the extents to form the buffer.
    sieve_size : int, default 9
        Maximum size of connected land pixels to sieve out from
        both low and high tide arrays to reduce noise over water.
        Set to 0 to apply no sieving.

    Returns
    -------
    coastal_mask : xarray.DataArray
        Boolean mask of the continuous coastal buffer zone.
    """
    # Create water frequency map (high = mostly wet)
    water_mask = mndwi.where(~inland_water_mask) > 0
    freq = water_mask.where(mndwi.notnull()).mean(dim="time")

    # Define extents with water represented as True
    high_tide = freq > 0.1  # includes pixels that are rarely wet
    low_tide = freq >= 0.9  # includes pixels that are mostly wet

    # Sieve water arrays to remove isolated wet pixels
    if sieve_size > 0:
        high_tide_cleaned = xr.apply_ufunc(
            sieve,
            high_tide.astype("int16"),
            kwargs={"size": sieve_size, "connectivity": 8},
            keep_attrs=True,
        ).astype(bool)
    
        low_tide_cleaned = xr.apply_ufunc(
            sieve,
            low_tide.astype("int16"),
            kwargs={"size": sieve_size, "connectivity": 8},
            keep_attrs=True,
        ).astype(bool)
    else:
        high_tide_cleaned = high_tide
        low_tide_cleaned = low_tide

    # Dilate the extents
    dilated_high_tide = mask_cleanup(
        high_tide_cleaned, mask_filters=[("dilation", dilation_pixels)]
    )
    dilated_low_tide = mask_cleanup(
        ~low_tide_cleaned, mask_filters=[("dilation", dilation_pixels)]
    )

    # Generate final buffer via intersection
    coastal_mask = dilated_high_tide & dilated_low_tide

    return coastal_mask, freq