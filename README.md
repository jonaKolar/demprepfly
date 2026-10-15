# demprepfly

Functions to prepare digital elevation models (DEM) for UAV mission planning on challenging terrain glaciers. The software loads and saves a DEM as geotif. Optional processing steps, in the order they are applied, are:

- **Reproject** the DEM horizontally (e.g. from lat / lon to a metric CRS) and / or vertically (e.g. between geoid and ellipsoid heights).
- **Fill holes** (e.g. crevasses or moulins) marked by a polygon shapefile, so the UAV doesn't fly down into them.
- **Detect crevasses automatically** using a morphological closing, and fill them the same way.
- **Fill nodata pixels** of the input DEM by interpolation.
- **Smooth** the DEM with a Gaussian filter of a given radius.
- **Vertically shift** the DEM by a constant offset.

Each step is optional and only runs if the relevant argument is provided, so the tool can be used for a single step (e.g. just filling holes) or several at once.

Only the first band of the input GeoTIFF is used. Nodata pixels are ignored during smoothing and are written back as nodata in the output (unless `--fill-nodata` is used). The output keeps the data type and nodata value of the input. Optionally, a mask of all pixels that got filled can be saved alongside the output.

All distances (`--detect-maxsize`, `--detect-enlarge`, `--fill-max-search-distance`, `--smoothing-radius`) are given in coordinate units of the DEM after reprojection, usually meters. As the reprojection runs first, a DEM in a lat / lon CRS can be reprojected to a metric CRS before the other steps. A warning is shown if filling, crevasse detection or smoothing run on a non metric CRS.

Can be imported and used in a script (see `usage.ipynb`), or run from the command line.

# Command line usage

```
python demprepfly.py input.tif output.tif [options]
```

**Required arguments**

| Argument | Description |
|---|---|
| `input_tif` | Path to the input DEM GeoTIFF. |
| `output_tif` | Path to save the output DEM GeoTIFF. |

**Output**

| Argument | Description | Default |
|---|---|---|
| `--save-mask` | Also save the mask of filled pixels (1 = filled, 0 = unchanged) as a GeoTIFF next to `output_tif`, with `_mask` added before the extension (e.g. `output_mask.tif`). | off |

**Reprojection**

Runs before all other steps. The horizontal CRS of the input is taken from `input_tif`. Under the hood, PROJ is used (via pyproj). If PROJ doesn't have the required geoid grid installed, it is downloaded (internet access needed on first use). If the grid can't be found, an error is raised instead of silently skipping the vertical shift.

| Argument | Description | Default |
|---|---|---|
| `--to-horiz-epsg` | EPSG code of the target horizontal CRS, e.g. `2056` (swiss LV95). | `None` (no horizontal reprojection) |
| `--from-vert-epsg` | Height reference of the input: EPSG code of a vertical CRS, e.g. `3855` (EGM2008) or `5729` (swiss LHN95), EPSG code of a non vertical CRS for heights above its ellipsoid, e.g. `4979` (WGS84), or `ellipsoid` for heights above the ellipsoid of the input CRS. Must be used together with `--to-vert-epsg`. | `None` (no vertical transformation) |
| `--to-vert-epsg` | Target height reference, same options as `--from-vert-epsg`. `ellipsoid` refers to the ellipsoid of the target CRS. Must be used together with `--from-vert-epsg`. | `None` (no vertical transformation) |

> **Ellipsoidal heights.** There are no vertical EPSG codes for ellipsoidal heights. They can be selected in two ways:
> - The EPSG code of a non vertical (geographic or projected) CRS means heights above the ellipsoid of *that* CRS, independent of the grid. E.g. `4979` (or `4326`) gives WGS84 ellipsoidal heights, also on a swiss LV95 grid.
> - `ellipsoid` means heights above the ellipsoid of the corresponding horizontal CRS (the target CRS for `--to-vert-epsg`, i.e. `--to-horiz-epsg`, or the input CRS if not given). This is not always WGS84: e.g. UTM32N (EPSG:32632) uses the WGS84 ellipsoid, but the swiss LV95 (EPSG:2056) uses Bessel 1841. In the Alps, heights above the two ellipsoids differ by about 50 m.
> - Examples: 
WGS84: EPSG:4326 or EPSG:4979 
ETRS89: EPSG:4258, 
CH1903+ / LV95: EPSG:4150 or EPSG:2056
CH1903 / LV03: EPSG:21781

>**Geoidal heights.** Examples:
EGM96 Geoid: EPSG:5773
EGM2008 Geoid: EPSG:3855
LN02 geoid: EPSG:5728
LHN95 geoid: EPSG:5729

**Hole filling from mask**

| Argument | Description | Default |
|---|---|---|
| `--fill-mask` | Path to polygon shapefile marking areas to remove and fill by interpolation. If omitted, no areas from a shapefile are filled. | `None` |

**Hole filling from crevasse detection**

Crevasses are detected by applying a morphological closing with a circular kernel, which fills depressions narrower than the kernel. Pixels where the original DEM lies more than the minimum depth below the closed DEM are flagged as crevasses, enlarged slightly to cover their edges, and then filled like the areas from `--fill-mask`. Detection runs only if both `--detect-maxsize` and `--detect-mindepth` are given; giving only one of them is an error.

| Argument | Description | Default |
|---|---|---|
| `--detect-maxsize` | Maximum crevasse width to detect, in coordinate units (usually meters). | `0` (no detection) |
| `--detect-mindepth` | Minimum crevasse depth to detect, in the units of `input_tif`. | `0` (no detection) |
| `--detect-enlarge` | Radius (in coordinate units) by which detected crevasses are enlarged before filling. | 10% of `detect-maxsize / 2` (rounded up to whole pixels) |

**Hole filling parameters**

| Argument | Description | Default |
|---|---|---|
| `--fill-nodata` | Also fill nodata pixels of the input DEM by interpolation. | off |
| `--fill-method` | Interpolation method for filling: `idw` (inverse distance weighting), `laplace` (smoothest surface matching the edge heights) or `biharmonic` (smoothest surface matching edge heights and slopes, continues the surrounding slope). See [`prep_dem`](#prep_dem) for details. | `laplace` |
| `--fill-max-search-distance` | Only for `--fill-method idw`: max distance (in coordinate units) to search in all directions for values to interpolate from. | `100` |
| `--fill-smoothing-iterations` | Only for `--fill-method idw`: number of 3x3 smoothing passes run on the filled section. | `3` |

**Smoothing**

| Argument | Description | Default |
|---|---|---|
| `--smoothing-radius` | Standard deviation (in coordinate units) of the Gaussian smoothing of the whole DEM. Nodata pixels are ignored. `0` disables smoothing. | `0.0` |

**Vertical shift**

| Argument | Description | Default |
|---|---|---|
| `--vshift` | Constant offset added to all values, in the units of `input_tif`. Use a negative value to lower the DEM. | `0.0` |

Use `python demprepfly.py -h` for the full help text.

### Examples

Fill crevasses from a shapefile and smooth the DEM with a 2 m Gaussian:

```
python demprepfly.py input.tif output.tif --fill-mask crevasses.shp --smoothing-radius 2.0
```

Automatically detect and fill depressions/crevasses up to 10 m wide and at least 0.5 m deep and save the detected areas as mask:

```
python demprepfly.py input.tif output.tif --detect-maxsize 10 --detect-mindepth 0.5 --save-mask
```

Fill nodata holes in the DEM:

```
python demprepfly.py input.tif output.tif --fill-nodata
```

Convert a DEM with WGS84 ellipsoidal heights in UTM 32N to EGM96 geoid heights and lower it by 20 m (e.g., to account for surface lowering through ice melt):

```
python demprepfly.py input.tif output.tif --from-vert-epsg 4326 --to-vert-epsg 5773 --vshift -20
```

Reproject a lat / lon DEM to swiss LV95 before detecting crevasses, so that the distances are in meters:

```
python demprepfly.py input.tif output.tif --to-horiz-epsg 2056 --detect-maxsize 10 --detect-mindepth 0.5
```

Reproject a DEM with swiss LHN95 heights to a swiss LV95 grid with WGS84 ellipsoidal heights (e.g. for a UAV flying on GNSS heights):

```
python demprepfly.py input.tif output.tif --to-horiz-epsg 2056 --from-vert-epsg 5729 --to-vert-epsg 4979
```


# Usage in python scripts
`demprepfly.py` must be in the same directory as your script. Then you can import functions:
```python
from demprepfly import prep_dem, reproj_dem, detect_crevasses, inpaint_dem
```

#### `prep_dem`

```python
prep_dem(input_tif, output_tif=None, save_mask=False, fill_mask=None, fill_nodata=False,
         detect_maxsize=0, detect_mindepth=0, detect_enlarge=None, fill_method="laplace",
         fill_max_search_distance=100, fill_smoothing_iterations=3, smoothing_radius=0.0, vshift=0.0,
         to_horiz_epsg=None, from_vert_epsg=None, to_vert_epsg=None)
```

Applies multiple preparation steps to a DEM GeoTIFF. Optional steps, in the order they are applied:

1. **Reprojection**: horizontal and / or vertical (see [`reproj_dem`](#reproj_dem))
2. **Fill holes**: remove and interpolate areas loaded from a shapefile
3. **Detect crevasses**: automatically detect crevasses, then remove and interpolate them
4. **Smoothing**: Gaussian smoothing
5. **Vertical shift**: constant offset

Distances (`detect_maxsize`, `detect_enlarge`, `fill_max_search_distance`, `smoothing_radius`) are given in coordinate units of the DEM after reprojection (usually meters) and converted to pixels internally. As the reprojection runs first, this allows to change away from a lat / lon CRS before working on the data. A warning is issued if filling, crevasse detection or smoothing run on a non metric CRS.

| Argument | Description | Default |
|---|---|---|
| `input_tif` | Path to DEM GeoTIFF. Only the first band is used. | |
| `output_tif` | Path to save the DEM GeoTIFF. If `None`, no file is saved. | `None` |
| `save_mask` | If `True`, the final mask is saved as uint8 GeoTIFF (1 = filled) next to `output_tif`, with `_mask` added before the extension (e.g. `dem.tif` → `dem_mask.tif`). Requires `output_tif`. | `False` |
| `fill_mask` | Path to a polygon shapefile marking areas to be removed and filled by interpolation. | `None` |
| `fill_nodata` | If `True`, nodata (and NaN) pixels of `input_tif` are also filled by interpolation. | `False` |
| `detect_maxsize` | Maximum width of crevasses to detect, in coordinate units. Detection only runs if both `detect_maxsize` and `detect_mindepth` are > 0. | `0` |
| `detect_mindepth` | Minimum depth of crevasses to detect, in the units of `input_tif`. | `0` |
| `detect_enlarge` | Radius by which detected crevasses are enlarged before filling, in coordinate units. If `None`, 10% of `detect_maxsize / 2` is used. | `None` |
| `fill_method` | Interpolation method for filling the masked areas, see below. | `"laplace"` |
| `fill_max_search_distance` | Only for `fill_method="idw"`: maximum distance to search in all directions for values to interpolate from, in coordinate units. | `100` |
| `fill_smoothing_iterations` | Only for `fill_method="idw"`: number of 3x3 smoothing passes run on the filled section. | `3` |
| `smoothing_radius` | Standard deviation of the Gaussian filter applied after filling, in coordinate units. | `0.0` |
| `vshift` | Offset added to all values, in the units of `input_tif`. Use negative values to lower the DEM. | `0.0` |
| `to_horiz_epsg` | EPSG code of the target horizontal CRS (e.g. `4326`). The source CRS is taken from `input_tif`. `None` means no horizontal reprojection. | `None` |
| `from_vert_epsg` | Height reference of `input_tif`: EPSG code of a vertical CRS (e.g. `3855` for EGM2008 or `5729` for swiss LHN95), EPSG code of a non vertical CRS for heights above its ellipsoid (e.g. `4979` for WGS84), or `"ellipsoid"` for heights above the ellipsoid of the CRS of `input_tif`. | `None` |
| `to_vert_epsg` | Target height reference, same options as `from_vert_epsg`. `"ellipsoid"` refers to the ellipsoid of the target horizontal CRS, which is not always WGS84 (e.g. Bessel 1841 for the swiss EPSG:2056). Use `4979` for WGS84 ellipsoidal heights on any grid. If both vertical codes are `None`, no vertical transformation is done. | `None` |

Fill methods:

| `fill_method` | Description |
|---|---|
| `"idw"` | Inverse distance weighting from the hole edges, followed by smoothing passes (rasterio `fillnodata`). Fast, but tends to give flat "tents" across wider holes on slopes. Holes wider than `fill_max_search_distance` are not completely filled. |
| `"laplace"` | Smoothest surface matching the heights at the hole edges (membrane). No overshoot, but does not continue the slope of the surroundings. |
| `"biharmonic"` | Smoothest surface matching heights and slopes at the hole edges (thin plate, skimage `inpaint_biharmonic`). Continues the surrounding slope, may overshoot at noisy edges. Slow for very large holes. |

Laplace and biharmonic fill all masked areas completely.

**Returns** a tuple `(dem, mask)`:
- `dem`: 2D `xarray.DataArray` of the prepared DEM (float32, NaN for nodata) with georeferencing (in the target CRS if reprojected). It is written to `output_tif` in the dtype and with the nodata value of `input_tif`.
- `mask`: 2D boolean `xarray.DataArray`, `True` where pixels were removed and filled by interpolation (from `fill_mask`, detected crevasses and, if `fill_nodata`, nodata pixels).

**Raises**
- `ValueError` if `save_mask` is `True` but `output_tif` is not provided.
- `ValueError` if only one of `detect_maxsize` and `detect_mindepth` is > 0.
- `ValueError` if `fill_method` is not one of `"idw"`, `"laplace"` or `"biharmonic"`.
- `ValueError` if only one of `from_vert_epsg` and `to_vert_epsg` is given.
- `ValueError` if a reprojection is requested but `input_tif` has no CRS.
- `pyproj.exceptions.ProjError` if the vertical transformation is not available, e.g. because the required geoid grid can't be found or downloaded.

#### `reproj_dem`

```python
reproj_dem(dem, from_horiz_epsg=None, to_horiz_epsg=None, from_vert_epsg=None, to_vert_epsg=None,
           n_samples=10, resampling=Resampling.bilinear)
```

Reprojects a DEM horizontally and vertically.

The vertical transformation is computed with pyproj on a coarse grid of `n_samples` x `n_samples` points of the DEM (spanning its full extent), linearly interpolated to all pixels and added to the DEM. This is accurate as long as the vertical shift varies smoothly (e.g. geoid models). The horizontal reprojection is then done with `rio.reproject`.

There are no vertical EPSG codes for ellipsoidal heights. They can be selected in two ways:
- `"ellipsoid"`: heights above the ellipsoid of the corresponding horizontal CRS (source CRS for `from_vert_epsg`, target CRS for `to_vert_epsg`). E.g. WGS84 for EPSG:32632, but Bessel 1841 for the swiss EPSG:2056 (heights differ by ~50 m in the Alps).
- EPSG code of a non vertical (geographic or projected) CRS: heights above the ellipsoid of that CRS, independent of the horizontal CRS. E.g. `4979` (or `4326`) for WGS84 ellipsoidal heights, also on a grid in EPSG:2056.

| Argument | Description | Default |
|---|---|---|
| `dem` | 2D rioxarray DataArray of elevations, with x / y coordinates in `from_horiz_epsg`. | |
| `from_horiz_epsg` | EPSG code of the horizontal CRS of `dem` (e.g. `32632`). If given, it overwrites the CRS of `dem`. `None` uses the CRS of `dem`. | `None` |
| `to_horiz_epsg` | EPSG code of the target horizontal CRS (e.g. `4326`). `None` means no horizontal reprojection. | `None` |
| `from_vert_epsg` | Height reference of `dem`: EPSG code of a vertical CRS (e.g. `3855` for EGM2008), EPSG code of a non vertical CRS for heights above its ellipsoid (e.g. `4979` for WGS84), or `"ellipsoid"` for heights above the ellipsoid of the source horizontal CRS. | `None` |
| `to_vert_epsg` | Target height reference, same options as `from_vert_epsg`. `"ellipsoid"` refers to the ellipsoid of the target horizontal CRS (`to_horiz_epsg`, or the source CRS if not given). If both vertical codes are `None`, no vertical transformation is done. | `None` |
| `n_samples` | Number of sample points per axis for computing the vertical shift. | `10` |
| `resampling` | Resampling method for the horizontal reprojection. | `Resampling.bilinear` |

**Returns** the reprojected DEM as rioxarray DataArray.

**Raises**
- `ValueError` if only one of `from_vert_epsg` and `to_vert_epsg` is given.
- `ValueError` if a transformation is requested but `dem` has no CRS and `from_horiz_epsg` is not given.
- `pyproj.exceptions.ProjError` if the vertical transformation is not available, e.g. because the required geoid grid can't be found or downloaded.

#### `detect_crevasses`

```python
detect_crevasses(dem, radius, enlarge=None, depth=1.0)
```

Detects crevasses (local depressions) in a DEM using a morphological closing.

A closing with a circular kernel fills depressions smaller than the radius. Pixels where the original DEM lies more than `depth` below the closed DEM are flagged as crevasses, and the resulting mask is dilated (enlarged) to also cover crevasse edges.

| Argument | Description | Default |
|---|---|---|
| `dem` | 2D array of elevation values. | |
| `radius` | Radius of the circular closing kernel, in pixels. Depressions up to about `2 * radius` wide are detected. | |
| `enlarge` | Radius of the circular kernel used to dilate the detected mask, in pixels. `None` uses 10% of `radius`. | `None` |
| `depth` | Minimum depth for a depression to count as a crevasse, in the units of the `dem` values. | `1.0` |

**Returns** a boolean array with the same shape as `dem`, `True` where crevasses were detected.

#### `inpaint_dem`

```python
inpaint_dem(dem, mask, method="laplace")
```

Fills masked pixels of a DEM by Laplace or biharmonic inpainting.

NaN pixels outside the mask (nodata) are temporarily set to the value of the nearest valid pixel, so that holes bordering nodata areas can be filled, and are set back to NaN afterwards.

| Argument | Description | Default |
|---|---|---|
| `dem` | 2D array of elevation values. | |
| `mask` | 2D boolean array, `True` where pixels are filled. | |
| `method` | `"laplace"` for the smoothest surface matching the heights at the hole edges (membrane), `"biharmonic"` for the smoothest surface matching heights and slopes (thin plate). | `"laplace"` |

**Returns** a 2D float array with the masked pixels filled.

See `usage.ipynb` for a full walkthrough.

## Requirements

- `rasterio`
- `rioxarray`
- `xarray`
- `pyproj` (>= 3.4)
- `geopandas`
- `numpy`
- `scipy`
- `scikit-image`
