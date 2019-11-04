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

## Special instructions for `Spyder`:
Change the backend to automatic:

```
Tools > preferences > IPython console > Graphics > Graphics backend > Backend: Automatic
```
(from [StackOverflow](https://stackoverflow.com/questions/23585126/how-do-i-get-interactive-plots-again-in-spyder-ipython-matplotlib))

Then close and open Spyder.

## Sampling raster data
Review [Sampling raster data](http://www.qgistutorials.com/en/docs/3/sampling_raster_data.html) to see our objective.

### Download and extract the data:

- [us.tmax_nohads_ll_20190501_float.tif](http://www.qgistutorials.com/downloads/us.tmax_nohads_ll_20190501_float.tif)
- [2018_Gaz_ua_national.zip](http://www.qgistutorials.com/downloads/2018_Gaz_ua_national.zip)
- [tl_2018_us_county.zip](http://www.qgistutorials.com/downloads/tl_2018_us_county.zip)

Open `Spyder` and create a new document in this repo named `sample_raster.py`.

### Import libraries
We are going to again use `geopandas` and  `descartes` in this lab, but also `rasterio`, `rasterstats`, `matplotlib`, and others. To start with:
```
import rasterio
```

### Load the data in python
We are going to reuse this string that gives the path to the tmax file you downloaded so let's store it in a variable.
When you use a string more than once, it is better to create a variable to store the value because it is very easy
for the value to become out of sync if, say, you want to change it and only change one of the instances.
```
tmax_path='/Users/aaryno/classes/gist604b/fall-2019-online/data/us.tmax_nohads_ll_20190501_float.tif'
```
Now, load the tif using `rasterio`:
```
dataset = rasterio.open(tmax_path)
```
Do some data exploration:
```
dataset.width
dataset.height
dataset.count
dataset.bounds
dataset.transform
dataset.transform * (dataset.width, dataset.height)
dataset.crs
dataset.indexes
band1 = dataset.read(1)
band1
```
Plot the raster:
```
from rasterio.plot import show
show(dataset)
```

### read csv
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
gaz = pandas.read_csv('/Users/aaryno/classes/gist604b/fall-2019-online/data/2018_Gaz_ua_national.txt',delimiter='\t')
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
Next, specify a coordinate reference system as a `dict`, since geopandas wants a dict like this when we construct our
spatial data frame:
```
crs = {'init': 'epsg:4326'}
```
Finally, create the spatial dataframe using the table, the geometries, and the crs:
```
gaz_geo = geopandas.GeoDataFrame(gaz, crs=crs, geometry=geometry)
```
Explore:
```
type(gaz_geo)
gaz_geo.head()
gaz_geo.plot()
```

### Plot them together:
```
rasterio.plot.show((dataset, 1))
ax = mpl.pyplot.gca()

patches = [PolygonPatch(feature) for feature in gaz_df]
ax.add_collection(mpl.collections.PatchCollection(patches))

mpl.figure(1)
import matplotlib.pyplot as plt
from matplotlib import interactive
interactive(True)
plt.figure(1)
gaz_geo.plot()
```

from matplotlib import pyplot
fig, ax = pyplot.subplots(1, figsize=(12, 12))
show(dataset, ax=ax)
gaz_geo.plot(ax=ax)
pyplot.show()

dataset.bounds.left

b = dataset.bound
dataset_env = geopandas.GeoSeries(shapely.geometry.Polygon([(dataset.bounds.left, dataset.bounds.bottom), (dataset.bounds.right, dataset.bounds.bottom), (dataset.bounds.right, dataset.bounds.top), (dataset.bounds.left, dataset.bounds.top),(dataset.bounds.left, dataset.bounds.bottom) ]))

df1 = geopandas.GeoDataFrame({'geometry': dataset_env, 'a':[1]})
dataset_env.plot()

gaz_48 = geopandas.sjoin(gaz_geo, df1)

gaz_48.plot()



# Sample
xy = [xy for xy in zip(gaz_48['INTPTLONG'], gaz_48.INTPTLAT)]

samples = [z for z in dataset.sample(xy=xy)]


gaz_48['tmax'] = samples


import rasterstats
counties_shapefile='/Users/aaryno/classes/gist604b/fall-2019-online/data/tl_2018_us_county/tl_2018_us_county.shp'
counties = geopandas.read_file(counties_shapefile)
counties.plot()

county_stats = rasterstats.zonal_stats(counties_shapefile, tmax_path, stats="count mean", all_touched=True)

type(county_stats)


tmax_mean = [x['mean'] for x in county_stats]
tmax_count = [x['count'] for x in county_stats]

counties['tmax_mean'] = tmax_mean


counties.to_file('/Users/aaryno/classes/gist604b/fall-2019-online/data/tl_2018_us_county/tl_2018_us_county_temps.shp')

counties.plot(column='tmax_mean',cmap='OrRd', scheme='quantiles')

from matplotlib import pyplot
fig, ax = pyplot.subplots(1, figsize=(12, 12))

ax.set_xlim([dataset.bounds.left, dataset.bounds.right])
ax.set_ylim([dataset.bounds.bottom, dataset.bounds.top])
counties.plot(column='tmax_mean', ax=ax)
gaz_geo.plot(ax=ax,marker='+')
pyplot.show()
