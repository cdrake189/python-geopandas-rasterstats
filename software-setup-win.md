# Python Geo Software Setup (Windows)

The power and flexibility of open source software means you can develop amazing new features or software using foundational
open source libraries, including libraries that handle projections, file format conversion, etc. As these libraries mature,
new features and versions come out but dependent software may be "pinned" with dependencies to old versions. Unfortunately,
one of the challenges in working with open source software is setting up a proper environment in which libraries are all
compatible with each other. This short tutorial shows you how to set up a working Python 3 environment in Windows 2010. 
Versions are listed below:

## Install OSGeo4w
- [Install OSGeo4W](https://trac.osgeo.org/osgeo4w/) - Choose 64-bit version for `network install` as long as your OS supports it
- Assuming OSGeo4W is installed in `c:\osgeo4w64`, add the following to your environment:
  - Add to the `Path`: `c:\osgeo4w64\bin` (confirm directory exists and contains libraries; e.g., `gdal111.dll`)
  - Add a new environment variable for rasterio: `GDAL_DATA` with the value: `c:\osgeo4w64\share\epsg_csv` (confirm directory exists)

## Install Anaconda for Python 3.7
- [Download Anaconda] (https://www.anaconda.com/distribution/)

## Conda environment setup
In your start menu, look for `Anaconda prompt` and open it. This will open a terminal screen with a prompt like this:
```
(base) c:\Users\aaryn> 
```
The `base` represents your conda _environment_. For this class we are going to create a _new_ environment and install our libraries there. Note that this means special care needs to be taken to make other Anaconda-installed software work with our environment which won't be covered in this class. Therefore, for future Anaconda-related projects (e.g., `Sypder` and `Jupyter`) we will launch them from the `Anaconda prompt` rather than through the Windows Start Menu.

### Add the `conda-forge` channel
`conda` is a package manager, which means you can type `conda install X` and it will install the library `X` as long as it 
resides in a repository (or _channel_) that conda knows about. Since `Anaconda` is a selected package of program and 
libraries optimized for data science, there is a special _channel_ for `conda`-compatible libraries. There are some
repositories we want to use that have not been `conda`-blessed but, whenever possible, we want `conda` to obtain
libraries directly from the `conda-forge` channel before looking at other repositories. At your anaconda prompt type:
```
conda config --add channels conda-forge
```

### Create a new environment
Set up a new conda environment with:
```
conda create --name python37 python=3.7
```
This creates a new python environment that uses python 3.7 (that's what the `python=3.7` means). It will be named `python37`. We will need to activate it in order to use it:
```
conda activate python37
```

At this point you have two environments installed. You have the default, or `base` environment and you have the `python37` environment. Type `conda env list` to see your list of installed environments.

If you close your shell and re-open it, if your conda environment is not active, you can activate it with 
`conda activate python37`. Rememeber that in the `Anaconda prompt` the environment will be listed in your prompt 
(i.e., `(base)` or `(python37)`. Of course, the name of the environment can be anything (i.e., you can have an 
environment named `bob` that used `python=3.7`).

For more details on conda environments: [Mananging conda environments](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)

### Install libraries
Now that we have created the environment we can install our libraries. In this class we will need `geopandas`, `descartes`, `shapely`, `rasterio`, and `raster_stats`. All but `raster_stats` can be installed from conda, although `rasterio` needs special handling.

First, install the easy ones:
```
conda install -y shapely geopandas descartes
```

Next, install `rasterio`, which has a special dependency. [`rasterio`](https://rasterio.readthedocs.io/) allows us to read and write rasters.
```
conda install -u rasterio proj4=5.2.0
```
This "pins" the `proj4` library at version `5.2.0`. Without that pin it would install the latest version which, at the time 
of this lab (November, 2019), is not compatible with `rasterio`.

Now, install `rasterstats`, which is not in the `conda-forge` channel ([rasterstats doc](https://pythonhosted.org/rasterstats/). `rasterstats` allows Zonal Statistics (e.g., summarize raster pixels based on polygons)
```
pip install rasterstats
```

Finally, install `spyder` and `jupyter`, which are GUIs that we will use to run python and edit our code interactively. 
`Anaconda` installed these but, as mentioned previously, it is not straightforward to align them to use our `conda env`
so we will be running these directly from the command line (i.e., `Anaconda prompt`)

```
conda install spyder jupyter
```

