#!/usr/bin/env python3
"""Command line tool to remove and interpolate crevasse areas from a GeoTIFF DEM.

Usage:
    python fill_crevasses.py input.tif input.shp output.tif [options]
"""

import argparse
import numpy as np


def fill_dem(input_tif, input_shp, output_tif=None, smoothing_radius=0.0,
                    offset=0.0, filling_max_search_distance=1000, smoothing_iterations=3):
    """Removes and interpolates areas from geotif indicated by shapefile.
    Afterwards, a Gaussian smoothening and a vertical offset is applied.

    Args:
        input_tif: Path to geotif.
        input_shp: Path to crevasse mask. Polygon shapefile marking areas to be removed and filled.
        output_tif: Path to save geotif. If None, no file is saved.
        smoothing_radius: Radius for gaussian filter after filling. In pixels.
        offset: Offset added to all values. Use negative values to lower the DEM. In units of the input.tif content.
        filling_max_search_distance: The maximum number of pixels to search in all directions to find values to interpolate from.
        smoothing_iterations: The number of 3x3 smoothing filter passes to run on the filled section.
    """

    from geopandas import read_file as gpd_read_file
    from rasterio import open as rasterio_open
    from rasterio.features import rasterize
    from rasterio.fill import fillnodata
    from scipy.ndimage import gaussian_filter

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

    # read crevasse mask
    gdf = gpd_read_file(input_shp).to_crs(crs)  # read crevasse shape file and transform to the raster CRS
    mask = rasterize([(geom, 1) for geom in gdf.geometry], out_shape=(height, width),
                      transform=transform, fill=0)[np.newaxis, :, :]
    mask = np.array(mask, dtype=np.bool)

    # fill masked values
    #data = np.where(mask, np.nan, data)
    data = fillnodata(data, mask=~mask, max_search_distance=filling_max_search_distance, smoothing_iterations=smoothing_iterations)

    # smoothening
    data = gaussian_filter(data, sigma=smoothing_radius)

    # apply offset
    data = data + offset

    data = data.astype(profile["dtype"])

    # Save with identical georeferencing and metadata
    if output_tif is not None:
        with rasterio_open(output_tif, "w", **profile) as dst:
            dst.write(data)

    return data


def parse_args():
    parser = argparse.ArgumentParser(
        description="Remove and interpolate crevasse areas from a GeoTIFF DEM, "
                     "using a polygon shapefile mask, then smooth and offset the result."
    )
    parser.add_argument("input_tif", help="Path to input GeoTIFF DEM.")
    parser.add_argument("input_shp", help="Path to polygon shapefile marking areas to remove and fill.")
    parser.add_argument("output_tif", help="Path to save the output GeoTIFF.")
    parser.add_argument("--smoothing-radius", type=float, default=0.0,
                         help="Radius (in pixels) for Gaussian smoothing after filling. Default: 0.0")
    parser.add_argument("--offset", type=float, default=0.0,
                         help="Vertical offset added to all values (negative to lower the DEM). Default: 0.0")
    parser.add_argument("--max-search-distance", type=float, default=1000,
                         help="Maximum search distance (in pixels) for gap filling. Default: 100.0")
    parser.add_argument("--smoothing_iterations", type=float, default=3,
                             help="The number of 3x3 smoothing filter passes to run on the filled section.")
    return parser.parse_args()


def main():
    args = parse_args()
    fill_dem(
        input_tif=args.input_tif,
        input_shp=args.input_shp,
        output_tif=args.output_tif,
        smoothing_radius=args.smoothing_radius,
        offset=args.offset,
        filling_max_search_distance=args.max_search_distance,
    )
    print(f"Saved DEM to {args.output_tif}")


if __name__ == "__main__":
    main()