# Python Spatial Programming with Geopandas, Rasterio, and Rasterstats
## Assignment

### Background
We previously worked with python within the framework of QGIS but python is most flexible when used in stand-alone mode, which allows us to decouple the QGIS framework libraries from the base libraries of python and explicitly add libraries of our own. As with other languages, the power of Python comes from the plethora of high-quality libraries that extend the base functionality of the core python programming language. In the realm of spatial python there are many libraries but we will start with a couple: `geopandas`, `shapely`, `rasterio`, and `rasterstats`.

Reference:
- [Shapely docs](https://shapely.readthedocs.io/en/stable/manual.html)
- [Geopanda docs](http://geopandas.org/)
- [Rasterio docs](https://pythonhosted.org/rasterstats/manual.html)
- [Rasterstats docs](https://pythonhosted.org/rasterstats/manual.html)

The objective of this lab is to reproduce one of the QGIS Tutorials you did previously:
- [Sampling raster data](http://www.qgistutorials.com/en/docs/3/sampling_raster_data.html)

## Deliverables
- `sample_raster.py`
- `sample_raster.png`
- `zonal_stats_county_temps.png`

## Install any new libraries
You will need to install at least `rasterstats` for this assignment. Add it to your `geo37` environment.

## Sampling raster data
Review [Sampling raster data](http://www.qgistutorials.com/en/docs/3/sampling_raster_data.html) to see our objective.

### Download and extract the data to a single directory:
Note that the paths in this assignment assume that the tl_2018_us_county shapefile is in the same directory as the other data and is _not_ in a subdirectory.

- [us.tmax_nohads_ll_20190501_float.tif](https://www.qgistutorials.com/downloads/us.tmax_nohads_ll_20190501_float.tif)
- [2018_Gaz_ua_national.zip](https://www.qgistutorials.com/downloads/2018_Gaz_ua_national.zip)
- [tl_2018_us_county.zip](https://www.qgistutorials.com/downloads/tl_2018_us_county.zip)

### Open `Spyder` and import libraries

We are going to again use `geopandas` and  `descartes` in this lab, but also `rasterio`, `rasterstats`, `matplotlib`, and others. To start with:
```
import os
import rasterio
```

### Load the data in python
We are going to reuse this string that gives the path to the tmax file you downloaded so let's store it in a variable.
When you use a string more than once, it is better to create a variable to store the value because it is very easy
for the value to become out of sync if, say, you want to change it and only change one of the instances.

For this assignment we assume all the data is in the same directory, which we'll label `base_dir`. Update this as needed for your environment.
```
base_dir = '/Users/aaryn/gist604b/python/rasterstats'
```
No open the raster tmax:
```
tmax_path=os.path.join(base_dir,'us.tmax_nohads_ll_20190501_float.tif')
```
Now, load the tif using `rasterio`:
```
tmax = rasterio.open(tmax_path)
```
Do some data exploration:
```
tmax.width
tmax.height
tmax.count
tmax.bounds
tmax.transform
tmax.transform * (tmax.width, tmax.height)
tmax.crs
tmax.indexes
band1 = tmax.read(1)
band1
```
Plot the raster:
```
from rasterio.plot import show
show(tmax)
```

### Read csv
Import libraries:
```
import pandas 
import shapely
import geopandas
import matplotlib as mpl
from descartes import PolygonPatch
```
Read CSV using pandas and specifying a tab (`\t`) delimiter:
```
gaz = pandas.read_csv(os.path.join(base_dir,'2018_Gaz_ua_national.txt'),delimiter='\t')
```
It does a poor job with the column names (one of them has a newline in it), so let's fix it:
```
gaz.columns = gaz.columns.str.strip()
```
Now explore:
```
gaz.head()
```
### Create points from csv:
This is a python way of creating an array from a loop. This creates an array of `shapely` geometries based on the two
columns, `INTPTLONG` and `INTPTLAT`:
```
geometry = [shapely.geometry.Point(xy) for xy in zip(gaz['INTPTLONG'], gaz.INTPTLAT)]
```
Explore:
```
type(geometry)
geometry[0]
```
Finally, create the spatial dataframe using the table, the geometries, and the crs:
```
gaz_geo = geopandas.GeoDataFrame(gaz, geometry=geometry, crs='epsg:4326')
```
Explore:
```
type(gaz_geo)
gaz_geo.head()
gaz_geo.plot()
```

### Plot them together:
When we create a plot we can get access to the plot object, including its axes, and manipulate them. Additionally, we can 
use the same plot and add data to it. Note that for this to work in spyder you will need to set your `iPython` ->  `Graphics` -> `Backend` to `Automatic`
```
rasterio.plot.show((tmax, 1))
ax = mpl.pyplot.gca()

gaz_geo.plot(ax=ax)
```
Note the extents change when the `gaz_geo` dataset is added to the plot. Let's zoom back in by changing the axis limits. Fix the axes to keep the extent snapped to the boundaries of the `tmax` raster
```
ax.set_xlim([tmax.bounds.left, tmax.bounds.right])
ax.set_ylim([tmax.bounds.bottom, tmax.bounds.top])
```
### Clip the points to the extent of the raster
Before we sample the points we will need to clip the points to the extent of the raster. Otherwise you will encounter errors
in the sampling. To clip points to the raster we need first create a polygon boundary of the raster and use a spatial
join to intersect the points with the polygon. Unfortunately, creating a polygon from the extent of the raster does not
seem to be very straight forward. We will first get the extents of the raster, then manually construct a geometry. In order
to use the `geopandas` `sjoin` we will then need to create a geopandas SpatialDataFrame. This section details those steps.

The `rasterio` raster has a `bounds` attribute which gives the boundaries:
```
tmax.bounds
```
To create a `shapely.geometry`:
```
shapely.geometry.Polygon([(left, bottom), (right, bottom), (right, top), (left, top), (left, bottom)])
```
which, subbing in the attributes from the `tmax.bounds`:
```
poly = shapely.geometry.Polygon([(tmax.bounds.left, tmax.bounds.bottom), (tmax.bounds.right, tmax.bounds.bottom), (tmax.bounds.right, tmax.bounds.top), (tmax.bounds.left, tmax.bounds.top),(tmax.bounds.left, tmax.bounds.bottom)])
poly
```
Geopandas wants the `Polygon` in a `GeoSeries`:
```
tmax_env_series = geopandas.GeoSeries(poly)
```
Finally, to create the `GeoDataFrame` we need a `geometry` and a `DataFrame`
```
tmax_env_gdf = geopandas.GeoDataFrame({'geometry': tmax_env_series, 'a':[1]})
tmax_env_gdf
tmax_env_gdf.head()
tmax_env_gdf.plot()
```
Now we can intersect the points with the polygons. This is doing using the geopandas `sjoin` method. See 
[geopandas doc](http://geopandas.org/mergingdata.html#spatial-joins) for more information.
```
gaz_48 = geopandas.sjoin(gaz_geo, tmax_env_gdf)
gaz_48
gaz_48.shape
gaz_48.plot()
```
If you did it right, then there should only be points within the outline of the `tmax` raster (i.e., in the continental US).

### Sample the points
To sample using the `rasterio` library we will need to create a sequence of `x,y` pairs. See [rasterio docs](https://rasterio.readthedocs.io/en/latest/api/rasterio._io.html#rasterio._io.DatasetReaderBase.sample).
```
xy = [xy for xy in zip(gaz_48['INTPTLONG'], gaz_48.INTPTLAT)]
```
Next, sample the raster with the `x,y` pairs:
```
tmax.sample(xy=xy)
```
You'll see it's a `generator object`. To convert it to a list (which we want so we can add a new column to our points 
GeoDataFrame, `gaz_48`), just call the `list()` function on it:
```
samples = list(tmax.sample(xy=xy))
```
`samples` should be the same length as the number of records in `gaz_48` and in the same order. 
```
len(samples)
gaz_48.shape
```
But the type of each sample is a `numpy.ndarray` - We want a simpler number. 
```
type(samples[0])
```
We need to unnest:
```
gaz_48['tmax'] = [item for unnested in samples for item in unnested]
```
Now we can append the `samples` list to the `gas_48` DataFrame:
```
gaz_48['tmax'] = samples
```
## Zonal stats
In the QGIS Tutorial, the tutorial switches gears and performs zonal stats on of the `mean` max temperature within each
county. This is a relatively straightforward task with `rasterstats` library, which has a `zonal_stats` function. The `zonal_stats` takes filenames as arguments. Continue to add commands to your `sample_raster.py` file.

First, import the counties shapefile:
```
import rasterstats

counties_path=os.path.join(base_dir,'tl_2018_us_county.shp')
counties = geopandas.read_file(counties_path)
counties.plot()
```
That's hard to see since we are interested in the continental US
```
ax = mpl.pyplot.gca()
counties.plot(ax=ax)
ax.set_xlim([tmax.bounds.left, tmax.bounds.right])
ax.set_ylim([tmax.bounds.bottom, tmax.bounds.top])
```
Take a minute to think about the resolutions and scales of these two datasets, the grid `tmax` and the vector, `counties`. If you look at the size of `tl_2018_us_county.shp` on disk you will see it is around 127.3 MB. On the other hand, the `tmax` tif is only 41KB. 
Additionally, the two data sources are in different projections. Let's reproject the data and do some vector generalization. When we reproject to EPSG:4326 (lat/long, WGS84), we'll use the distance 0.01 degrees as the simplication tolerance.
```
counties_4326 = counties.to_crs({'init': 'epsg:4326'}).simplify(0.001)
counties_4326_path = os.path.join(base_dir,'tl_2018_us_county_4326.shp')
counties_4326.to_file(counties_4326_path)
```
Next, perform the zonal stats. We need to set `all_touched` = `True` (see [rasterstats docs](https://pythonhosted.org/rasterstats/manual.html#rasterization-strategy) for discussion.
```
county_stats = rasterstats.zonal_stats(counties_4326_path, tmax_path, stats="mean", all_touched=True, nodata=-999)
```
Dig in and see what it gives us:
```
type(county_stats)
county_stats[0]
type(county_stats[0])
```
This gives us a `list` of `dict`s with the `mean` for each county. This is almost exactly what we want. But we also
want it joined to the counties GeoDataFrame but it needs to be converted into a list of numbers. 

```
tmax_mean = [x['mean'] for x in county_stats]
type(tmax_mean)
tmax_mean[0]
len(tmax_mean)
```
This looks like what we want; just like with the samples, we will create a new column on `counties` to store the 
mean `tmax`:
```
counties['tmax_mean'] = tmax_mean
```
Next, save the file:
```
counties.to_file(os.path.join(base_dir,'tl_2018_us_county_temps.shp'))
```

### Plot the counties based on the mean max temp
Plot the counties using the `tmax_mean` column:
```
counties.plot(column='tmax_mean')
```
That's not easy to see so let's change the axes:
```
ax.set_xlim([tmax.bounds.left, tmax.bounds.right])
ax.set_ylim([tmax.bounds.bottom, tmax.bounds.top])
counties.plot(column='tmax_mean', ax=ax)
```
Let's stretch the color map:
```
counties.plot(column='tmax_mean',cmap='OrRd', scheme='quantiles')
```
To finish, let's make a plot containing both the `gaz_48` points and the counties symbolized by mean max temperature.
```
from matplotlib import pyplot
fig, ax = pyplot.subplots(1, figsize=(12, 12))

ax.set_xlim([tmax.bounds.left, tmax.bounds.right])
ax.set_ylim([tmax.bounds.bottom, tmax.bounds.top])
counties.plot(column='tmax_mean', ax=ax)
gaz_48.plot(ax=ax,marker='+',color='red')
pyplot.show()
```

Save this as `sample_raster.png`

## Deliverable
Look at the [geopandas docs](http://geopandas.org/reference.html#geopandas.GeoDataFrame.plot) to finish up. Recreate the last map containing the counties symbolized based on `tmax_mean` and the `gaz_48` urban areas but change:
1) the color scheme of the counties
2) the color of the gaz urban areas
3) the symbol of the gaz urban areas

Create a screenshot of your final map and name it `zonal_stats_county_temps.png `

[![DOI](https://zenodo.org/badge/219385180.svg)](https://zenodo.org/badge/latestdoi/219385180)
