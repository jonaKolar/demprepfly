#!/usr/bin/env python3
"""Command line tool to prepare a DEM GeoTIFF: fill holes, smooth, vertically shift,
and/or convert between geoid and ellipsoid heights.

Usage:
    python prepdem.py input.tif output.tif [options]
"""

import argparse

import numpy as np
from geopandas import read_file as gpd_read_file
from rasterio import open as rasterio_open
from rasterio.warp import reproject, Resampling
from rasterio.features import rasterize
from rasterio.fill import fillnodata
from scipy.ndimage import gaussian_filter


def prepdem(input_tif, output_tif=None, fill_mask=None, fill_max_search_distance=1000,
            fill_smoothing_iterations=3, smoothing_radius=0.0, vshift=0.0,
            to_ellipsoid=False, to_geoid=False, geoid_grid=None):
    """Applies multiple preparation steps to a digital elevation model (DEM) tif file. Optional steps are:
        - fill holes: remove and interpolate areas loaded from shape file
        - smoothening: applies Gaussian smoothening
        - vertical shift: constant offset
        - to geoid / ellipsoid: applies vertical shift defined in loaded grid file

    Args:
        input_tif: Path to DEM geotif.
        output_tif: Path to save DEM geotif. If None, no file is saved. Defaults to None.
        fill_mask: Path to .shp file mask. Polygon shapefile marking areas to be removed and filled by interpolation. Defaults to None.
        fill_max_search_distance: The maximum number of pixels to search in all directions to find values to interpolate from. Defaults to 1000.
        fill_smoothing_iterations: The number of 3x3 smoothing filter passes to run on the filled section. Defaults to 3.
        smoothing_radius: Radius for gaussian filter after filling. In pixels. Defaults to 0.0.
        vshift: Offset added to all values. Use negative values to lower the DEM. In units of the input.tif content. Defaults to 0.0.
        to_ellipsoid: If True geoid height is added to the data. geoid_grid must be provided. Defaults to False.
        to_geoid: If True geoid height is subtracted from the data. geoid_grid must be provided. Defaults to False.
        geoid_grid: Path to gridded geoid height as .agr or .tif file. Defaults to None.

    Raises:
        ValueError: If to_ellipsoid or to_geoid is True but geoid_grid is not provided.

    Returns:
        data
    """

    # Read dem data
    with rasterio_open(input_tif) as src:
        data = src.read()  # shape: (bands, rows, cols)
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
                          transform=transform, fill=0)[np.newaxis, :, :]
        mask = np.array(mask, dtype=bool)

        # fill masked values
        data = fillnodata(data, mask=~mask, max_search_distance=fill_max_search_distance,
                           smoothing_iterations=fill_smoothing_iterations)

    # smoothening
    if smoothing_radius > 0.0:
        data = gaussian_filter(data, sigma=smoothing_radius)

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
                source=src.read(),
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

    data = data.astype(profile["dtype"])

    # Save with identical georeferencing and metadata
    if output_tif is not None:
        with rasterio_open(output_tif, "w", **profile) as dst:
            dst.write(data)

    return data


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare a DEM GeoTIFF: optionally fill masked areas, apply Gaussian "
                     "smoothing, a vertical shift, and/or geoid/ellipsoid conversion."
    )
    parser.add_argument("input_tif", help="Path to input DEM GeoTIFF.")
    parser.add_argument("output_tif", help="Path to save the output DEM GeoTIFF.")

    parser.add_argument("--fill-mask", default=None,
                         help="Path to polygon shapefile marking areas to remove and fill. Default: None (no filling).")
    parser.add_argument("--fill-max-search-distance", type=float, default=1000,
                         help="Max pixel search distance for gap filling. Default: 1000")
    parser.add_argument("--fill-smoothing-iterations", type=int, default=3,
                         help="Number of 3x3 smoothing passes on the filled section. Default: 3")

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

    prepdem(
        input_tif=args.input_tif,
        output_tif=args.output_tif,
        fill_mask=args.fill_mask,
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