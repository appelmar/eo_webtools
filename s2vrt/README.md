# s2vrt
Generate GDAL VRT datasets for Sentinel 2 data

## Description
This commmand line tool generates GDAL VRT datasets for Sentinel 2 (L1C) imagery and optionally converts JPEG200 files to cloud optimized GeoTIFFs and builds overviews.
For fast interactive mapping, tiled GeoTIFF files with overview levels 2, 4, 8, 16, 32 should be prefered even though they might consume more disk space.


## Usage
```
Input arguments:
  -h [ --help ]         print usage
  -o [ --overview ] arg overview levels which will be generated for input 
                        images
  --cog                 convert input Sentinel 2 subdatasets to cloud-optimized
                        GeoTIFFs 
  -v [ --verbose ]      verbose output
  --input-file arg
  --output-file arg
```

# Examples

**Build a simple VRT dataset to index all Sentinel 2 bands in one image**
```
s2vrt  S2A_MSIL1C_20170815T102021_N0205_R065_T32TMR_20170815T102513.SAFE/MTD_MSIL1C.xml test.vrt
```


**Convert a S2 scene to cloud-optimized GeoTIFF files including overviews**  
```
s2vrt -v -o "2 4 8 16 32" --cog S2A_MSIL1C_20170815T102021_N0205_R065_T32TMR_20170815T102513.SAFE/MTD_MSIL1C.xml S2A_MSIL1C_20170815T102021_N0205_R065_T32TMR_20170815T102513.vrt
```

# Build instructions
(see Dockerfile)
