# demprepfly

Functions to prepare DEMs for UAV mission planning.

## `demprepfly.py`

Prepares a DEM GeoTIFF for flight planning. It can:

- **Fill holes** (e.g. crevasses or moulins) marked by a polygon shapefile, so the UAV doesn't fly down into them.
- **Smooth** the DEM with a Gaussian filter of a given radius.
- **Vertically shift** the DEM by a constant offset.
- **Convert between geoid and ellipsoid heights** using a gridded geoid file.

Each step is optional and only runs if the relevant argument is provided, so the tool can be used for a single step (e.g. just filling holes) or several at once.

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

**Hole filling**

| Argument | Description | Default |
|---|---|---|
| `--fill-mask` | Path to polygon shapefile marking areas to remove and fill by interpolation. If omitted, no filling is done. | `None` |
| `--fill-max-search-distance` | Max number of pixels to search in all directions for values to interpolate from. | `1000` |
| `--fill-smoothing-iterations` | Number of 3x3 smoothing passes run on the filled section. | `3` |

**Smoothing**

| Argument | Description | Default |
|---|---|---|
| `--smoothing-radius` | Radius (in pixels) for Gaussian smoothing of the whole DEM. `0` disables smoothing. | `0.0` |

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

Just convert an ellipsoid-height DEM to geoid heights:

```
python demprepfly.py input.tif output.tif --to-geoid --geoid-grid geoid.tif
```

Fill holes only, with a larger search distance (e.g. for large holes (measured in pixels):

```
python demprepfly.py input.tif output.tif --fill-mask crevasses.shp --fill-max-search-distance 2000
```

### Python usage

```python
from prepdem import prepdem

data = prepdem(
    input_tif="input.tif",
    output_tif="output.tif",
    fill_mask="crevasses.shp",
    smoothing_radius=5.0,
    vshift=-20.0,
)
```

See `usage.ipynb` for a full walkthrough.

## Requirements

- `rasterio`
- `geopandas`
- `numpy`
- `scipy`
