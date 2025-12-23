'''
Generates a clipped and projected cropland cover map based on the MODIS
MCD12Q1 Type 5 land-cover classification.
'''

import os
import numpy as np
import h5py
import rasterio
import py4eos
import pyproj
from rasterio.warp import calculate_default_transform, reproject, Resampling

hdf = py4eos.read_file(snakemake.input[0], platform = 'MODIS')
hdf

lc_raster = hdf.to_rasterio('LC_Type5', filename = '', driver = 'MEM')

# NOTE: 6933 is the EPSG code for the global EASE-Grid 2.0 projection
new_transform, width, height = calculate_default_transform(
    lc_raster.crs, pyproj.CRS(6933), *lc_raster.shape, *lc_raster.bounds)

# See more at: https://rasterio.readthedocs.io/en/stable/topics/resampling.html

output_raster = rasterio.open(
    snakemake.output[0], 'w+', count = 1, width = width, height = height,
    dtype = np.uint8, crs = pyproj.CRS(6933), transform = new_transform)

# Writing the new array to the file
reproject(
    source = rasterio.band(lc_raster, 1),
    destination = rasterio.band(output_raster, 1),
    resampling = Resampling.mode,
    src_nodata = 255,
    dst_nodata = 255)

lc_map_1km = output_raster.read(1)

# Create a binary array where Croplands==1, 0 for everything else
croplands = np.where(lc_map_1km == 7, 1, 0)

output_raster.write(croplands, 1)
output_raster.close()
