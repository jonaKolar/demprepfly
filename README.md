# flightDEM
Functions to prepare DEMs for UAV mission planning

```fill_dem```: Function to fill holes like crevasses or moulins in DEMs. Additionally, the DEM can be shifted up or down by a specified amount and a the result can be Gaussian-smoothed by a given radius. This function can be imported and used in a script (see ```usage.ipynb```) or it can be used from the command line. Command line usage:

~~~
python fill_dem.py input.tif input.shp output.tif --smoothing-radius 5.0 --offset -20.0
~~~

Use ```python fill_crevasses.py -h``` for help
