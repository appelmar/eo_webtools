# eo_reduce
A command line tool to reduce Earth observation image collections over time.

## Description

The eo_reduce command line tool takes a collection of input images as input and applies a reducer function on pixel time series.
It warps individual scenes to the specified projection, spatial window, and resolution and afterwards reduces all the images.

All images must have the same number of bands and the same band data types but may have different spatial footprints. Mosaicing is performed 
automatically based on the given output properties.  




### Usage
```
eo_reduce OPTIONS SOURCE_IMAGE_1 [SOURCE_IMAGE_2]*
eo_reduce OPTIONS SOURCE_DIR
```

**Options:**
```
  -h [ --help ]              print usage
  -o [ --output-file ] arg   path of output file
  -b [ --bands ] arg         Bands of the input dataset, defaults to all 
                             available bands
  --t_win arg                coordinates of the target window expressed in the 
                             target reference system  with order (xmin, xmax, 
                             ymin, ymax).
  --t_size arg               size of the target image in pixels with order 
                             (width, height).
  --t_srs arg                target reference system as GDAL readable string.
  -f [ --func ] arg          reducer function, currently supported are 'mean, 
                             'min', 'max', 'median' 
  -m [ --memory ] arg (=256) memory in MB that is available for caching
  -v [ --verbose ]           verbose output
```


### Examples
**Derive mean NDVI image from Landsat 7 imagery**
```
eo_reduce -f "mean" --t_win="799500 809500 789000 799000" --t_size="334 334" --t_srs="EPSG:32636" -o test_reduce_landsat.tif ../data/Landsat/
```



## Build Instructions
(see Dockerfile)


### Dependencies
* [gdal](http://www.gdal.org/)
* [proj](http://proj4.org/)
* [geos](https://trac.osgeo.org/geos/)
* [Boost](http://www.boost.org/) `program_options` and `filesystem` libraries
* [cmake](https://cmake.org/) >= 3.5
* A C++11 compiler
