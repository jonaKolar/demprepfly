#!/usr/bin/env python3
"""Command line tool to prepare a DEM GeoTIFF: fill holes, smooth, vertically shift,
and/or convert between geoid and ellipsoid heights.

Usage:
    python prepdem.py input.tif output.tif [options]
"""

import argparse
import os

import numpy as np
from geopandas import read_file as gpd_read_file
from rasterio import open as rasterio_open
from rasterio.warp import reproject, Resampling
from rasterio.features import rasterize
from rasterio.fill import fillnodata
from scipy.ndimage import gaussian_filter
from skimage.morphology import closing, disk
from scipy.ndimage import gaussian_filter, binary_dilation


def prepdem(input_tif, output_tif=None, save_mask=False, fill_mask=None, fill_nodata=False, detect_maxsize=0, detect_mindepth=0, detect_enlarge=None, fill_max_search_distance=1000,
            fill_smoothing_iterations=3, smoothing_radius=0.0, vshift=0.0,
            to_ellipsoid=False, to_geoid=False, geoid_grid=None):
    """Applies multiple preparation steps to a digital elevation model (DEM) tif file. Optional steps are:
        - fill holes: remove and interpolate areas loaded from shape file
        - detect crevasses: automatically detect crevasses, then remove and interpolate them
        - smoothening: applies Gaussian smoothening
        - vertical shift: constant offset
        - to geoid / ellipsoid: applies vertical shift defined in loaded grid file

    Args:
        input_tif: Path to DEM geotif. Only the first band is used.
        output_tif: Path to save DEM geotif. If None, no file is saved. Defaults to None.
        save_mask: If True, the final mask is saved as uint8 geotif (1 = filled) next to output_tif, with _mask
            added before the extension (e.g. dem.tif -> dem_mask.tif). Requires output_tif. Defaults to False.
        fill_mask: Path to .shp file mask. Polygon shapefile marking areas to be removed and filled by interpolation. Defaults to None.
        fill_nodata: If True, nodata (and NaN) pixels of input_tif are also filled by interpolation. Defaults to False.
        detect_maxsize: Maximum width of crevasses to detect. In pixels. Detection only runs if both detect_maxsize and detect_mindepth are > 0. Defaults to 0.
        detect_mindepth: Minimum depth of crevasses to detect. In units of the input.tif content. Defaults to 0.
        detect_enlarge: Radius by which detected crevasses are enlarged before filling. In pixels. If None, 10% of detect_maxsize / 2 is used. Defaults to None.
        fill_max_search_distance: The maximum number of pixels to search in all directions to find values to interpolate from. Defaults to 1000.
        fill_smoothing_iterations: The number of 3x3 smoothing filter passes to run on the filled section. Defaults to 3.
        smoothing_radius: Radius for gaussian filter after filling. In pixels. Defaults to 0.0.
        vshift: Offset added to all values. Use negative values to lower the DEM. In units of the input.tif content. Defaults to 0.0.
        to_ellipsoid: If True geoid height is added to the data. geoid_grid must be provided. Defaults to False.
        to_geoid: If True geoid height is subtracted from the data. geoid_grid must be provided. Defaults to False.
        geoid_grid: Path to gridded geoid height as .agr or .tif file. Defaults to None.

    Raises:
        ValueError: If save_mask is True but output_tif is not provided.
        ValueError: If only one of detect_maxsize and detect_mindepth is > 0.
        ValueError: If to_ellipsoid or to_geoid is True but geoid_grid is not provided.

    Returns:
        tuple: (data, mask)
            data: 2D array of the prepared DEM, in the dtype of input_tif.
            mask: 2D boolean array, True where pixels were removed and filled by interpolation
                (from fill_mask, detected crevasses and, if fill_nodata, nodata pixels).
    """

    if save_mask and output_tif is None:
        raise ValueError("output_tif must be provided if save_mask is used.")

    if (detect_maxsize > 0) != (detect_mindepth > 0):
        raise ValueError("detect_maxsize and detect_mindepth must both be provided for crevasse detection.")

    # Read dem data
    with rasterio_open(input_tif) as src:
        data = np.array(src.read(1), dtype=np.float32)  # shape: (rows, cols), select 1st band
        transform = src.transform
        crs = src.crs
        height = src.height
        width = src.width
        profile = src.profile.copy()
        nodata = src.nodata

    if nodata is not None:
        data[data == nodata] = np.nan

    # load shape file mask, remove and interpolate
    if fill_mask is not None:
        gdf = gpd_read_file(fill_mask).to_crs(crs)  # read crevasse shape file and transform to the raster CRS
        mask = rasterize([(geom, 1) for geom in gdf.geometry], out_shape=(height, width),
                          transform=transform, fill=0)
        mask = np.array(mask, dtype=bool)
    else:
        mask = np.zeros(data.shape).astype(bool)

    if detect_maxsize > 0 and detect_mindepth > 0:
        mask_detected = detect_crevasses(data, radius=.5*detect_maxsize, depth=detect_mindepth, enlarge=detect_enlarge)
        mask = mask | mask_detected

    if fill_nodata:
        mask = mask | np.isnan(data)

    # fill masked values
    data = fillnodata(data, mask=~mask, max_search_distance=fill_max_search_distance,
                        smoothing_iterations=fill_smoothing_iterations)

    # smoothening (normalised convolution to ignore NaN)
    if smoothing_radius > 0.0:
        valid = ~np.isnan(data)
        data_smooth = gaussian_filter(np.where(valid, data, 0.0), sigma=smoothing_radius)
        weights = gaussian_filter(valid.astype(np.float32), sigma=smoothing_radius)
        with np.errstate(invalid="ignore", divide="ignore"):
            data = np.where(valid, data_smooth / weights, np.nan)


    # apply offset
    data += vshift

    # geoid / ellipsoid conversion
    if to_ellipsoid or to_geoid:
        if geoid_grid is None:
            raise ValueError("geoid_grid must be provided if to_ellipsoid or to_geoid is used.")

        with rasterio_open(geoid_grid) as src:
            geoid_data = np.zeros_like(data, dtype=np.float32)

            if src.crs is None:
                geoid_crs = 'EPSG:4326'
            else:
                geoid_crs = src.crs

            reproject(
                source=src.read(1),
                destination=geoid_data,
                src_transform=src.transform,
                src_crs=geoid_crs,
                dst_transform=transform,
                dst_crs=crs,
                resampling=Resampling.bilinear
            )

        if to_ellipsoid:
            data += geoid_data

        if to_geoid:
            data -= geoid_data

    if nodata is not None:
        data[np.isnan(data)] = nodata


    data = data.astype(profile["dtype"])

    # Save with identical georeferencing and metadata
    if output_tif is not None:
        with rasterio_open(output_tif, "w", **profile) as dst:
            dst.write(data[np.newaxis, :, :])

    if save_mask:
        root, ext = os.path.splitext(output_tif)
        mask_profile = profile.copy()
        mask_profile.update(dtype="uint8", count=1, nodata=None)
        mask_profile.pop("predictor", None)  # a floating point predictor is invalid for uint8
        with rasterio_open(f"{root}_mask{ext}", "w", **mask_profile) as dst:
            dst.write(mask.astype(np.uint8)[np.newaxis, :, :])

    return data, mask


def detect_crevasses(dem, radius, enlarge=None, depth=1.0):
    """Detects crevasses (local depressions) in a DEM using a morphological closing.

    A closing with a circular kernel fills depressions smaller than the radius.
    Pixels where the original DEM lies more than `depth` below the closed DEM are flagged as
    crevasses, and the resulting mask is dilated (enlarged) to also cover crevasse edges.

    Args:
        dem: 2D array of elevation values.
        radius: Radius of the circular closing kernel. In pixels. Depressions up to about
            2 * radius wide are detected.
        enlarge: Radius of the circular kernel used to dilate the detected mask. In pixels.
            Defaults to 10% * radius.
        depth: Minimum depth for a depression to count as a crevasse. In units of the dem
            values. Defaults to 1.0.

    Returns:
        Boolean array with the same shape as dem, True where crevasses were detected.
    """
    if enlarge is None:
        enlarge = np.ceil(.1 * radius)

    # fill crevasses using closing operation with circular kernel
    dem_closed = closing(dem, footprint=disk(radius))

    # detect crevasses by subtracting original
    mask = (dem_closed - dem) > depth

    # enlarge True regions in mask (detected crevasses)
    mask = binary_dilation(mask, structure=disk(enlarge))

    return mask



def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare a DEM GeoTIFF: optionally fill masked areas, apply Gaussian "
                     "smoothing, a vertical shift, and/or geoid/ellipsoid conversion."
    )
    parser.add_argument("input_tif", help="Path to input DEM GeoTIFF.")
    parser.add_argument("output_tif", help="Path to save the output DEM GeoTIFF.")
    parser.add_argument("--save-mask", action="store_true",
                         help="Also save the mask of filled pixels as <output_tif>_mask.tif.")

    parser.add_argument("--fill-mask", default=None,
                         help="Path to polygon shapefile marking areas to remove and fill. Default: None (no filling).")
    parser.add_argument("--fill-max-search-distance", type=float, default=1000,
                         help="Max pixel search distance for gap filling. Default: 1000")
    parser.add_argument("--fill-smoothing-iterations", type=int, default=3,
                         help="Number of 3x3 smoothing passes on the filled section. Default: 3")
    parser.add_argument("--fill-nodata", action="store_true",
                         help="Also fill nodata pixels of the input DEM by interpolation.")

    parser.add_argument("--detect-maxsize", type=float, default=0,
                         help="Maximum crevasse width (in pixels) for automatic detection. "
                              "Detection requires --detect-maxsize and --detect-mindepth > 0. Default: 0 (no detection)")
    parser.add_argument("--detect-mindepth", type=float, default=0,
                         help="Minimum crevasse depth (in DEM units) for automatic detection. Default: 0 (no detection)")
    parser.add_argument("--detect-enlarge", type=float, default=None,
                         help="Radius (in pixels) by which detected crevasses are enlarged. Default: 10%% of detect-maxsize / 2")

    parser.add_argument("--smoothing-radius", type=float, default=0.0,
                         help="Radius (in pixels) for Gaussian smoothing. 0 disables smoothing. Default: 0.0")

    parser.add_argument("--vshift", type=float, default=0.0,
                         help="Vertical offset added to all values (negative to lower the DEM). Default: 0.0")

    geoid_group = parser.add_mutually_exclusive_group()
    geoid_group.add_argument("--to-ellipsoid", action="store_true",
                              help="Add geoid height to the data (requires --geoid-grid).")
    geoid_group.add_argument("--to-geoid", action="store_true",
                              help="Subtract geoid height from the data (requires --geoid-grid).")
    parser.add_argument("--geoid-grid", default=None,
                         help="Path to gridded geoid height (.agr or .tif). Required with --to-ellipsoid/--to-geoid.")

    return parser.parse_args()


def main():
    args = parse_args()

    if (args.to_ellipsoid or args.to_geoid) and args.geoid_grid is None:
        raise SystemExit("Error: --geoid-grid is required when using --to-ellipsoid or --to-geoid.")

    if (args.detect_maxsize > 0) != (args.detect_mindepth > 0):
        raise SystemExit("Error: --detect-maxsize and --detect-mindepth must be used together.")

    prepdem(
        input_tif=args.input_tif,
        output_tif=args.output_tif,
        save_mask=args.save_mask,
        fill_mask=args.fill_mask,
        fill_nodata=args.fill_nodata,
        detect_maxsize=args.detect_maxsize,
        detect_mindepth=args.detect_mindepth,
        detect_enlarge=args.detect_enlarge,
        fill_max_search_distance=args.fill_max_search_distance,
        fill_smoothing_iterations=args.fill_smoothing_iterations,
        smoothing_radius=args.smoothing_radius,
        vshift=args.vshift,
        to_ellipsoid=args.to_ellipsoid,
        to_geoid=args.to_geoid,
        geoid_grid=args.geoid_grid,
    )
    print(f"Saved prepared DEM to {args.output_tif}")


if __name__ == "__main__":
    main()