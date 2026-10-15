#!/usr/bin/env python3
"""Command line tool to prepare a DEM GeoTIFF: reproject horizontally and / or vertically,
fill holes and detected crevasses, smooth and / or vertically shift.

Usage:
    python demprepfly.py input.tif output.tif [options]
"""

import argparse
import os
import warnings

import numpy as np
from rioxarray import open_rasterio
from xarray import DataArray
from geopandas import read_file as gpd_read_file
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.fill import fillnodata
from skimage.morphology import closing, disk
from skimage.restoration import inpaint_biharmonic
from scipy.ndimage import gaussian_filter, binary_dilation, distance_transform_edt
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve
from pyproj import CRS, Transformer, network
from pyproj.crs import CompoundCRS


FILL_METHODS = ("idw", "laplace", "biharmonic")


def prep_dem(input_tif, output_tif=None, save_mask=False, fill_mask=None, fill_nodata=False, detect_maxsize=0, detect_mindepth=0, detect_enlarge=None,
            fill_method="laplace", fill_max_search_distance=100, fill_smoothing_iterations=3, smoothing_radius=0.0, vshift=0.0, to_horiz_epsg=None, from_vert_epsg=None, to_vert_epsg=None):
    """Applies multiple preparation steps to a digital elevation model (DEM) tif file. Optional steps,
    in the order they are applied:
        - reprojection: horizontal and / or vertical reprojection (see reproj_dem)
        - fill holes: remove and interpolate areas loaded from shape file
        - detect crevasses: automatically detect crevasses, then remove and interpolate them
        - smoothening: applies Gaussian smoothening
        - vertical shift: constant offset

    Distances (detect_maxsize, detect_enlarge, fill_max_search_distance, smoothing_radius) are given in
    coordinate units of the DEM after reprojection (usually meters) and converted to pixels internally.
    As the reprojection runs first, this allows to change away from a lat / lon CRS before working on
    the data. A warning is issued if filling, crevasse detection or smoothing run on a non metric CRS.

    Args:
        input_tif: Path to DEM geotif. Only the first band is used.
        output_tif: Path to save DEM geotif. If None, no file is saved. Defaults to None.
        save_mask: If True, the final mask is saved as uint8 geotif (1 = filled) next to output_tif, with _mask
            added before the extension (e.g. dem.tif -> dem_mask.tif). Requires output_tif. Defaults to False.
        fill_mask: Path to .shp file mask. Polygon shapefile marking areas to be removed and filled by interpolation. Defaults to None.
        fill_nodata: If True, nodata (and NaN) pixels of input_tif are also filled by interpolation. Defaults to False.
        detect_maxsize: Maximum width of crevasses to detect. In coordinate units. Detection only runs if both detect_maxsize and detect_mindepth are > 0. Defaults to 0.
        detect_mindepth: Minimum depth of crevasses to detect. In units of the input.tif content. Defaults to 0.
        detect_enlarge: Radius by which detected crevasses are enlarged before filling. In coordinate units. If None, 10% of detect_maxsize / 2 is used. Defaults to None.
        fill_method: Interpolation method for filling the masked areas. Defaults to "laplace".
            - "idw": inverse distance weighting from the hole edges, followed by smoothing passes
              (rasterio fillnodata). Fast, but tends to give flat "tents" across wider holes on slopes.
              Holes wider than fill_max_search_distance are not completely filled.
            - "laplace": smoothest surface matching the heights at the hole edges (membrane). No overshoot,
              but does not continue the slope of the surroundings.
            - "biharmonic": smoothest surface matching heights and slopes at the hole edges (thin plate,
              skimage inpaint_biharmonic). Continues the surrounding slope, may overshoot at noisy edges.
              Slow for very large holes.
            Laplace and biharmonic fill all masked areas completely.
        fill_max_search_distance: Only for fill_method="idw". The maximum distance to search in all directions to find values to interpolate from. In coordinate units. Defaults to 100.
        fill_smoothing_iterations: Only for fill_method="idw". The number of 3x3 smoothing filter passes to run on the filled section. Defaults to 3.
        smoothing_radius: Standard deviation of the gaussian filter applied after filling. In coordinate units. Defaults to 0.0.
        vshift: Offset added to all values. Use negative values to lower the DEM. In units of the input.tif content. Defaults to 0.0.
        to_horiz_epsg: EPSG code of the target horizontal CRS (e.g. 4326). The source CRS is taken from input_tif.
            Defaults to None, which means no horizontal reprojection.
        from_vert_epsg: Height reference of input_tif: EPSG code of a vertical CRS (e.g. 3855 for EGM2008 heights
            or 5729 for swiss LHN95), EPSG code of a non vertical CRS for heights above its ellipsoid (e.g. 4979
            for WGS84), or "ellipsoid" for heights above the ellipsoid of the CRS of input_tif. Defaults to None.
        to_vert_epsg: Target height reference, same options as from_vert_epsg. "ellipsoid" refers to the ellipsoid
            of the target horizontal CRS (to_horiz_epsg, or the CRS of input_tif if not given), which is not
            always WGS84, e.g. Bessel 1841 for the swiss EPSG:2056. Use 4979 for WGS84 ellipsoidal heights
            on any grid (see reproj_dem).
            Defaults to None. If both vertical codes are None, no vertical transformation is done.

    Raises:
        ValueError: If save_mask is True but output_tif is not provided.
        ValueError: If only one of detect_maxsize and detect_mindepth is > 0.
        ValueError: If fill_method is not one of "idw", "laplace" or "biharmonic".
        ValueError: If only one of from_vert_epsg and to_vert_epsg is given.
        ValueError: If a reprojection is requested but input_tif has no CRS.
        pyproj.exceptions.ProjError: If the vertical transformation is not available, e.g. because
            the required geoid grid can't be found or downloaded.

    Returns:
        tuple: (dem, mask)
            dem: 2D xarray.DataArray of the prepared DEM (float32, NaN for nodata) with georeferencing
                (in the target CRS if reprojected).
                It is written to output_tif in the dtype and with the nodata value of input_tif.
            mask: 2D boolean xarray.DataArray, True where pixels were removed and filled by interpolation
                (from fill_mask, detected crevasses and, if fill_nodata, nodata pixels).
    """

    if save_mask and output_tif is None:
        raise ValueError("output_tif must be provided if save_mask is used.")

    if (detect_maxsize > 0) != (detect_mindepth > 0):
        raise ValueError("detect_maxsize and detect_mindepth must both be provided for crevasse detection.")

    if fill_method not in FILL_METHODS:
        raise ValueError(f"fill_method must be one of {FILL_METHODS}, got '{fill_method}'.")

    # Read dem data, select 1st band. masked=True sets nodata to NaN; the original dtype
    # and nodata value are kept in the encoding and restored by to_raster.
    dem = open_rasterio(input_tif, masked=True).sel(band=1)

    dem = reproj_dem(dem, to_horiz_epsg=to_horiz_epsg, from_vert_epsg=from_vert_epsg, to_vert_epsg=to_vert_epsg)

    # the distance based steps below assume a metric grid, e.g. not degrees of a lat / lon crs
    distance_steps = fill_mask is not None or fill_nodata or detect_maxsize > 0 or smoothing_radius > 0
    if distance_steps and dem.rio.crs is not None:
        units = {axis.unit_name for axis in CRS(dem.rio.crs).axis_info[:2]}
        if units != {"metre"}:
            warnings.warn(f"Filling, crevasse detection and smoothing run on a non metric CRS "
                          f"({dem.rio.crs}, units: {', '.join(sorted(units))}). Distance parameters are "
                          f"interpreted in these units and pixels may not be square on the ground. "
                          f"Consider reprojecting to a metric CRS with to_horiz_epsg.")

    # pixel size in coordinate units, to convert distances to pixels
    res_x, res_y = np.abs(dem.rio.resolution())
    res = np.mean([res_x, res_y])

    data = dem.values.astype(np.float32)

    # ---------- work on the data ----------

    # load shape file mask, remove and interpolate
    if fill_mask is not None:
        gdf = gpd_read_file(fill_mask).to_crs(dem.rio.crs)  # read crevasse shape file and transform to the raster CRS
        mask = rasterize([(geom, 1) for geom in gdf.geometry], out_shape=data.shape,
                          transform=dem.rio.transform(), fill=0)
        mask = np.array(mask, dtype=bool)
    else:
        mask = np.zeros(data.shape).astype(bool)

    if detect_maxsize > 0 and detect_mindepth > 0:
        radius_px = max(1, round(.5 * detect_maxsize / res))
        enlarge_px = None if detect_enlarge is None else round(detect_enlarge / res)
        mask_detected = detect_crevasses(data, radius=radius_px, depth=detect_mindepth, enlarge=enlarge_px)
        mask = mask | mask_detected

    if fill_nodata:
        mask = mask | np.isnan(data)

    # fill masked values
    if fill_method == "idw":
        data = fillnodata(data, mask=~mask, max_search_distance=fill_max_search_distance / res,
                            smoothing_iterations=fill_smoothing_iterations)
    else:
        data = inpaint_dem(data, mask, method=fill_method)

    # smoothening (normalised convolution to ignore NaN), sigma per axis (rows, cols) in pixels
    if smoothing_radius > 0.0:
        sigma = (smoothing_radius / res_y, smoothing_radius / res_x)
        valid = ~np.isnan(data)
        data_smooth = gaussian_filter(np.where(valid, data, 0.0), sigma=sigma)
        weights = gaussian_filter(valid.astype(np.float32), sigma=sigma)
        with np.errstate(invalid="ignore", divide="ignore"):
            data = np.where(valid, data_smooth / weights, np.nan)

    # apply offset
    data += vshift

    # ---------- place changed data in georeferenced xarrays as before ----------
    dem = dem.copy(data=data)
    mask = dem.copy(data=mask)
    mask.attrs, mask.encoding = {}, {}

    # Save with identical georeferencing, dtype and nodata value
    if output_tif is not None:
        dem.rio.to_raster(output_tif)

    if save_mask:
        root, ext = os.path.splitext(output_tif)
        mask.astype(np.uint8).rio.to_raster(f"{root}_mask{ext}")

    return dem, mask




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


def inpaint_dem(dem, mask, method="laplace"):
    """Fills masked pixels of a DEM by Laplace or biharmonic inpainting.

    NaN pixels outside the mask (nodata) are temporarily set to the value of the nearest valid pixel,
    so that holes bordering nodata areas can be filled, and are set back to NaN afterwards.

    Args:
        dem: 2D array of elevation values.
        mask: 2D boolean array, True where pixels are filled.
        method: "laplace" for the smoothest surface matching the heights at the hole edges (membrane),
            "biharmonic" for the smoothest surface matching heights and slopes (thin plate). Defaults to "laplace".

    Returns:
        2D float array with the masked pixels filled.
    """
    dem = np.array(dem, dtype=np.float64)
    nodata = np.isnan(dem) & ~mask
    known = ~np.isnan(dem) & ~mask
    if not mask.any() or not known.any():
        return dem.astype(np.float32)  # nothing to fill, or nothing to interpolate from

    # temporarily replace nodata outside the mask by the nearest valid value
    nearest = distance_transform_edt(~known, return_distances=False, return_indices=True)
    filled = dem[tuple(nearest)]
    filled[mask] = 0.0

    if method == "laplace":
        filled = _inpaint_laplace(filled, mask)
    elif method == "biharmonic":
        filled = inpaint_biharmonic(filled, mask, split_into_regions=True)
    else:
        raise ValueError(f"method must be 'laplace' or 'biharmonic', got '{method}'.")

    filled[nodata] = np.nan
    return filled.astype(np.float32)


def _inpaint_laplace(image, mask):
    """Solves the Laplace equation on the masked pixels, with the unmasked pixels as boundary values."""
    n_rows, n_cols = image.shape
    rows, cols = np.nonzero(mask)
    n = rows.size
    index = np.full(image.shape, -1)
    index[rows, cols] = np.arange(n)

    # 5 point stencil: number of neighbors * z_i - sum of masked neighbors = sum of known neighbors.
    # Neighbors outside the image are left out (zero slope at the image border).
    diag = np.zeros(n)
    rhs = np.zeros(n)
    a_rows, a_cols = [], []
    for d_row, d_col in ((-1, 0), (1, 0), (0, -1), (0, 1)):
        r, c = rows + d_row, cols + d_col
        inside = (r >= 0) & (r < n_rows) & (c >= 0) & (c < n_cols)
        diag += inside
        i, r, c = np.nonzero(inside)[0], r[inside], c[inside]
        is_masked = mask[r, c]
        a_rows.append(i[is_masked])
        a_cols.append(index[r[is_masked], c[is_masked]])
        np.add.at(rhs, i[~is_masked], image[r[~is_masked], c[~is_masked]])

    a_rows, a_cols = np.concatenate(a_rows), np.concatenate(a_cols)
    A = coo_matrix((-np.ones(a_rows.size), (a_rows, a_cols)), shape=(n, n)) + diags(diag)

    out = image.copy()
    out[rows, cols] = spsolve(A.tocsr(), rhs)
    return out





def reproj_dem(dem, from_horiz_epsg=None, to_horiz_epsg=None, from_vert_epsg=None, to_vert_epsg=None,
               n_samples=10, resampling=Resampling.bilinear):
    """Reprojects a DEM horizontally and vertically.

    The vertical transformation is computed with pyproj on a coarse grid of n_samples x n_samples
    points of the DEM (spanning its full extent), linearly interpolated to all pixels and added
    to the DEM. This is accurate as long as the vertical shift varies smoothly (e.g. geoid models).
    The horizontal reprojection is then done with rio.reproject.

    There are no vertical EPSG codes for ellipsoidal heights. They can be selected in two ways:
        - "ellipsoid": heights above the ellipsoid of the corresponding horizontal CRS (source CRS for
          from_vert_epsg, target CRS for to_vert_epsg). E.g. WGS84 for EPSG:32632, but Bessel 1841 for
          the swiss EPSG:2056 (heights differ by ~50 m in the Alps).
        - EPSG code of a non vertical (geographic or projected) CRS: heights above the ellipsoid of that
          CRS, independent of the horizontal CRS. E.g. 4979 (or 4326) for WGS84 ellipsoidal heights,
          also on a grid in EPSG:2056.

    Args:
        dem: 2D rioxarray DataArray of elevations, with x / y coordinates in from_horiz_epsg.
        from_horiz_epsg: EPSG code of the horizontal CRS of dem (e.g. 32632). If given, it overwrites
            the CRS of dem. Defaults to None, which uses the CRS of dem.
        to_horiz_epsg: EPSG code of the target horizontal CRS (e.g. 4326). Defaults to None, which
            means no horizontal reprojection.
        from_vert_epsg: Height reference of dem: EPSG code of a vertical CRS (e.g. 3855 for EGM2008
            heights), EPSG code of a non vertical CRS for heights above its ellipsoid (e.g. 4979 for WGS84),
            or "ellipsoid" for heights above the ellipsoid of the source horizontal CRS. Defaults to None.
        to_vert_epsg: Target height reference, same options as from_vert_epsg. "ellipsoid" refers to the
            ellipsoid of the target horizontal CRS (to_horiz_epsg, or the source CRS if not given).
            Defaults to None. If both vertical codes are None, no vertical transformation is done.
        n_samples: Number of sample points per axis for computing the vertical shift. Defaults to 10.
        resampling: Resampling method for the horizontal reprojection. Defaults to Resampling.bilinear.

    Raises:
        ValueError: If only one of from_vert_epsg and to_vert_epsg is given.
        ValueError: If a transformation is requested but dem has no CRS and from_horiz_epsg is not given.
        pyproj.exceptions.ProjError: If the vertical transformation is not available, e.g. because
            the required geoid grid can't be found or downloaded.

    Returns:
        Reprojected DEM as rioxarray DataArray.
    """
    if (from_vert_epsg is None) != (to_vert_epsg is None):
        raise ValueError("from_vert_epsg and to_vert_epsg must both be given or both be None. "
                         "Use 'ellipsoid' for ellipsoidal heights.")

    # source horizontal CRS: given or taken from dem. Target defaults to source (no reprojection).
    if from_horiz_epsg is not None:
        dem = dem.rio.write_crs(f"EPSG:{from_horiz_epsg}")
    elif dem.rio.crs is None:
        if to_horiz_epsg is not None or from_vert_epsg != to_vert_epsg:
            raise ValueError("dem has no CRS. Provide from_horiz_epsg.")
        return dem  # nothing to do

    from_horiz_crs = CRS(dem.rio.crs)
    to_horiz_crs = from_horiz_crs if to_horiz_epsg is None else CRS.from_epsg(to_horiz_epsg)

    # ------------------------ vertical transform ------------------------
    if from_vert_epsg != to_vert_epsg:
        def height_crs(horiz_crs, vert_epsg):
            """Returns the 3D CRS describing the heights, and the 2D CRS of its x / y coordinates."""
            if vert_epsg == "ellipsoid":  # heights above the ellipsoid of the horizontal CRS
                return horiz_crs.to_3d(), horiz_crs
            vert_crs = CRS.from_epsg(vert_epsg)
            if vert_crs.is_vertical:  # e.g. geoid heights
                return CompoundCRS(f"{horiz_crs.name} + {vert_crs.name}", [horiz_crs, vert_crs]), horiz_crs
            # non vertical CRS: heights above its ellipsoid, x / y in its geographic (lon / lat) CRS
            geodetic_crs = vert_crs.geodetic_crs
            return geodetic_crs.to_3d(), geodetic_crs

        from_crs_3d, from_xy_crs = height_crs(from_horiz_crs, from_vert_epsg)
        to_crs_3d, _ = height_crs(to_horiz_crs, to_vert_epsg)  # only heights are used from the output

        network.set_network_enabled(True)  # allows PROJ to download missing geoid grids
        # only_best: raise an error instead of silently falling back to a transformation without vertical shift
        transformer = Transformer.from_crs(from_crs_3d, to_crs_3d, always_xy=True, only_best=True)

        # sample grid spanning the full dem, heights taken from the dem
        xs = np.linspace(dem.x.values[0], dem.x.values[-1], n_samples)
        ys = np.linspace(dem.y.values[0], dem.y.values[-1], n_samples)
        z = dem.sel(x=xs, y=ys, method="nearest").values.astype(np.float64)
        z[np.isnan(z)] = np.nanmedian(z) if np.any(~np.isnan(z)) else 0.0
        x, y = np.meshgrid(xs, ys)

        # express sample positions in the x / y CRS of the source heights, if it differs from the dem
        if from_xy_crs != from_horiz_crs:
            x, y = Transformer.from_crs(from_horiz_crs, from_xy_crs, always_xy=True).transform(x, y)

        _, _, z_new = transformer.transform(x, y, z)
        dz = DataArray(z_new - z, coords={"y": ys, "x": xs}, dims=("y", "x"))

        # interpolate vertical shifts to all pixels and add them (copy keeps encoding / nodata)
        dz = dz.interp(x=dem.x, y=dem.y).values
        dem = dem.copy(data=(dem.values + dz).astype(dem.dtype))

    # ------------------------ horizontal transform ------------------------
    if to_horiz_crs != from_horiz_crs:
        dem = dem.rio.reproject(to_horiz_crs, resampling=resampling)

    return dem



def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare a DEM GeoTIFF: optionally reproject horizontally and/or vertically, fill masked "
                     "areas and detected crevasses, apply Gaussian smoothing and/or a vertical shift."
    )
    parser.add_argument("input_tif", help="Path to input DEM GeoTIFF.")
    parser.add_argument("output_tif", help="Path to save the output DEM GeoTIFF.")
    parser.add_argument("--save-mask", action="store_true",
                         help="Also save the mask of filled pixels as <output_tif>_mask.tif.")

    parser.add_argument("--fill-mask", default=None,
                         help="Path to polygon shapefile marking areas to remove and fill. Default: None (no filling).")
    parser.add_argument("--fill-method", choices=FILL_METHODS, default="laplace",
                         help="Interpolation method for filling: 'idw' (inverse distance weighting, rasterio fillnodata), "
                              "'laplace' (smoothest surface matching edge heights) or 'biharmonic' (smoothest surface "
                              "matching edge heights and slopes, continues the surrounding slope). Default: laplace")
    parser.add_argument("--fill-max-search-distance", type=float, default=100,
                         help="Only for --fill-method idw: max search distance (in coordinate units, usually meters) "
                              "for gap filling. Default: 100")
    parser.add_argument("--fill-smoothing-iterations", type=int, default=3,
                         help="Only for --fill-method idw: number of 3x3 smoothing passes on the filled section. Default: 3")
    parser.add_argument("--fill-nodata", action="store_true",
                         help="Also fill nodata pixels of the input DEM by interpolation.")

    parser.add_argument("--detect-maxsize", type=float, default=0,
                         help="Maximum crevasse width (in coordinate units, usually meters) for automatic detection. "
                              "Detection requires --detect-maxsize and --detect-mindepth > 0. Default: 0 (no detection)")
    parser.add_argument("--detect-mindepth", type=float, default=0,
                         help="Minimum crevasse depth (in DEM units) for automatic detection. Default: 0 (no detection)")
    parser.add_argument("--detect-enlarge", type=float, default=None,
                         help="Radius (in coordinate units) by which detected crevasses are enlarged. Default: 10%% of detect-maxsize / 2")

    parser.add_argument("--smoothing-radius", type=float, default=0.0,
                         help="Standard deviation (in coordinate units) of the Gaussian smoothing. 0 disables smoothing. Default: 0.0")

    parser.add_argument("--vshift", type=float, default=0.0,
                         help="Vertical offset added to all values (negative to lower the DEM). Default: 0.0")

    parser.add_argument("--to-horiz-epsg", type=int, default=None,
                         help="EPSG code of the target horizontal CRS, e.g. 2056. Runs before all other steps. "
                              "Default: None (no horizontal reprojection)")
    parser.add_argument("--from-vert-epsg", type=vert_epsg, default=None,
                         help="Height reference of the input: EPSG code of a vertical CRS, e.g. 3855 (EGM2008), "
                              "EPSG code of a non vertical CRS for heights above its ellipsoid, e.g. 4979 (WGS84), "
                              "or 'ellipsoid' for heights above the ellipsoid of the input CRS. Requires --to-vert-epsg.")
    parser.add_argument("--to-vert-epsg", type=vert_epsg, default=None,
                         help="Target height reference, same options as --from-vert-epsg. 'ellipsoid' refers to the "
                              "ellipsoid of the target CRS (--to-horiz-epsg, or the input CRS if not given), which is "
                              "not always WGS84, e.g. Bessel 1841 for EPSG:2056. Use 4979 for WGS84 ellipsoidal heights "
                              "on any grid. Requires --from-vert-epsg.")

    return parser.parse_args()


def vert_epsg(value):
    """argparse type for vertical CRS arguments: an EPSG code or 'ellipsoid'."""
    if value == "ellipsoid":
        return value
    try:
        return int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"must be an EPSG code or 'ellipsoid', got '{value}'")


def main():
    args = parse_args()

    if (args.from_vert_epsg is None) != (args.to_vert_epsg is None):
        raise SystemExit("Error: --from-vert-epsg and --to-vert-epsg must be used together.")

    if (args.detect_maxsize > 0) != (args.detect_mindepth > 0):
        raise SystemExit("Error: --detect-maxsize and --detect-mindepth must be used together.")

    prep_dem(
        input_tif=args.input_tif,
        output_tif=args.output_tif,
        save_mask=args.save_mask,
        fill_mask=args.fill_mask,
        fill_nodata=args.fill_nodata,
        detect_maxsize=args.detect_maxsize,
        detect_mindepth=args.detect_mindepth,
        detect_enlarge=args.detect_enlarge,
        fill_method=args.fill_method,
        fill_max_search_distance=args.fill_max_search_distance,
        fill_smoothing_iterations=args.fill_smoothing_iterations,
        smoothing_radius=args.smoothing_radius,
        vshift=args.vshift,
        to_horiz_epsg=args.to_horiz_epsg,
        from_vert_epsg=args.from_vert_epsg,
        to_vert_epsg=args.to_vert_epsg,
    )
    print(f"Saved prepared DEM to {args.output_tif}")


if __name__ == "__main__":
    main()