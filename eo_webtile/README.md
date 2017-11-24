# eo_webtile
Tiling of Earth observation images with on-the-fly band selection, color scaling, and reprojection using GDAL


## Description

The eo_webtile command line tool generates one or more tiles from a GDAL dataset. Its main use is to enable interactive web-mapping of Earth observation datasets such as Sentinel 2 data. The tool's functionality is similar to `gdal2tiles.py` but aims at generating tiles on-demand with changing visualization parameters such as color scaling and band selection. The tool performs the following steps:

1. Extract selected bands and scale the values to 8-bit unsigned integers
2. Reproject and resize the image based on the tile's resolution and spatial extent susing `gdalwarp`
3. Convert the image to PNG with alpha channel

To make the process as efficient as possible, it does not store intermediate files on disk, only step 3 creates the final PNG file. However, computation times for tile generation strongly depend on how the input dataset is organized. In general, we recommend storing input data as cloud-optimized (tiled) GeoTIFF files 
including overviews. 

## Getting Started

**Create a single RGB tile from of a multispectral VRT dataset**
```
eo_webtile  -b "1 2 3" -z 10 -x 578 -y 569  ~/Desktop/out.vrt  ~/Desktop
```

### Usage

`eo_webtile OPTION SOURCEFILE TARGETDIR`

**Options**:
```
-h [ --help ]              print usage
-z [ --zoom ] arg          Zoom level(s)
-x [ --tx ] arg            tile x cooordinate(s)
-y [ --ty ] arg            tile x cooordinate(s)
-b [ --bands ] arg         Bands of the input dataset, either one or three 
                            must be selected
--min arg                  minimum values that will be scaled to 0 for each 
                            band
--max arg                  maximum values that will be scaled to 255 for each
                            band
-m [ --memory ] arg        maximum memory in MiB for gdalwarp and GDAL cache                                (defaut is 256 MiB)
-v [ --verbose ]           verbose output
```


### Examples

**Generate a single tile for bands 4,3,1 from an input dataset**
``` 
eo_webtile  -b "4 3 1" -z 10 -x 578 -y 569  image.tif /tmp/tiles
```

**Generate six tiles with the same zoom level at the same time**
```
eo_webtile  -b "1 2 3" -z "13 13 13 13 13 13" 
            -x "4627 4628 4629 4627 4628 4629" 
            -y "4556 4556 4556 4557 4557 4557" 
             image.tif  /tmp/tiles
```





## Build Instructions
(see Dockerfile)


### Dependencies
* [gdal](http://www.gdal.org/)
* [proj](http://proj4.org/)
* [Boost](http://www.boost.org/) `program_options` and `filesystem` libraries
* [cmake](https://cmake.org/) >= 3.5
* A C++11 compiler