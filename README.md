# demprepfly

Functions to prepare DEMs for UAV mission planning.

## `demprepfly.py`

Prepares a DEM GeoTIFF for flight planning. It can:

- **Fill holes** (e.g. crevasses or moulins) marked by a polygon shapefile, so the UAV doesn't fly down into them.
- **Detect crevasses automatically** using a morphological closing, and fill them the same way.
- **Fill nodata pixels** of the input DEM by interpolation.
- **Smooth** the DEM with a Gaussian filter of a given radius.
- **Vertically shift** the DEM by a constant offset.
- **Convert between geoid and ellipsoid heights** using a gridded geoid file.

Each step is optional and only runs if the relevant argument is provided, so the tool can be used for a single step (e.g. just filling holes) or several at once.

Only the first band of the input GeoTIFF is used. Nodata pixels are ignored during smoothing and are written back as nodata in the output (unless `--fill-nodata` is used). Optionally, a mask of all filled pixels can be saved alongside the output.

Can be imported and used in a script (see `usage.ipynb`), or run from the command line.

### Command line usage

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

**Hole filling**

| Argument | Description | Default |
|---|---|---|
| `--fill-mask` | Path to polygon shapefile marking areas to remove and fill by interpolation. If omitted, no areas from a shapefile are filled. | `None` |
| `--fill-nodata` | Also fill nodata pixels of the input DEM by interpolation. | off |
| `--fill-max-search-distance` | Max number of pixels to search in all directions for values to interpolate from. | `1000` |
| `--fill-smoothing-iterations` | Number of 3x3 smoothing passes run on the filled section. | `3` |

**Crevasse detection**

Crevasses are detected by applying a morphological closing with a circular kernel, which fills depressions narrower than the kernel. Pixels where the original DEM lies more than the minimum depth below the closed DEM are flagged as crevasses, enlarged slightly to cover their edges, and then filled like the areas from `--fill-mask`. Detection runs only if both `--detect-maxsize` and `--detect-mindepth` are given; giving only one of them is an error.

| Argument | Description | Default |
|---|---|---|
| `--detect-maxsize` | Maximum crevasse width to detect, in pixels. | `0` (no detection) |
| `--detect-mindepth` | Minimum crevasse depth to detect, in the units of `input_tif`. | `0` (no detection) |
| `--detect-enlarge` | Radius (in pixels) by which detected crevasses are enlarged before filling. | 10% of `detect-maxsize / 2` (rounded up) |

**Smoothing**

| Argument | Description | Default |
|---|---|---|
| `--smoothing-radius` | Radius (in pixels) for Gaussian smoothing of the whole DEM. Nodata pixels are ignored. `0` disables smoothing. | `0.0` |

**Vertical shift**

| Argument | Description | Default |
|---|---|---|
| `--vshift` | Constant offset added to all values, in the units of `input_tif`. Use a negative value to lower the DEM. | `0.0` |

**Geoid / ellipsoid conversion**

| Argument | Description | Default |
|---|---|---|
| `--to-ellipsoid` | Add geoid height to the data (converts geoid heights to ellipsoid heights). Requires `--geoid-grid`. Mutually exclusive with `--to-geoid`. | off |
| `--to-geoid` | Subtract geoid height from the data (converts ellipsoid heights to geoid heights). Requires `--geoid-grid`. Mutually exclusive with `--to-ellipsoid`. | off |
| `--geoid-grid` | Path to a gridded geoid height file (`.agr` or `.tif`). Required if `--to-ellipsoid` or `--to-geoid` is used. | `None` |

Use `python demprepfly.py -h` for the full help text.

### Examples

Fill crevasses, smooth, and lower the DEM by 20 meters:

```
python demprepfly.py input.tif output.tif --fill-mask crevasses.shp --smoothing-radius 5.0 --vshift -20.0
```

Automatically detect crevasses up to 30 pixels wide and at least 2 m deep, and save the mask of filled pixels:

```
python demprepfly.py input.tif output.tif --detect-maxsize 30 --detect-mindepth 2.0 --save-mask
```

Fill nodata holes in the DEM:

```
python demprepfly.py input.tif output.tif --fill-nodata
```

Just convert an ellipsoid-height DEM to geoid heights:

```
python demprepfly.py input.tif output.tif --to-geoid --geoid-grid geoid.tif
```

Fill holes only, with a larger search distance (e.g. for large holes, measured in pixels):

```
python demprepfly.py input.tif output.tif --fill-mask crevasses.shp --fill-max-search-distance 2000
```

### Python usage

`prepdem()` returns a tuple `(data, mask)`: the prepared DEM as a 2D array (in the dtype of the input), and a 2D boolean array that is `True` where pixels were removed and filled.

```python
from demprepfly import prepdem

data, mask = prepdem(
    input_tif="input.tif",
    output_tif="output.tif",
    save_mask=True,
    fill_mask="crevasses.shp",
    detect_maxsize=30,
    detect_mindepth=2.0,
    smoothing_radius=5.0,
    vshift=-20.0,
)
```

`output_tif` can be omitted to only get the arrays back without saving a file (`save_mask` then can't be used).

See `usage.ipynb` for a full walkthrough.

## Requirements

- `rasterio`
- `geopandas`
- `numpy`
- `scipy`
- `scikit-image`
