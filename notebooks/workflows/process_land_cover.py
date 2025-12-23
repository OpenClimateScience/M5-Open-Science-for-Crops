'''
Generates a clipped and projected cropland cover map based on the MODIS
MCD12Q1 Type 5 land-cover classification.
'''

import glob
import os
import earthaccess
import numpy as np
import h5py
import rasterio
import py4eos
import pyproj
from rasterio.warp import calculate_default_transform, reproject, Resampling

auth = earthaccess.login()

TIME_PERIOD = ('2023-10-01', '2024-09-30')

# The bounding box for our study area
bbox = (1.5, 34.0, 8.0, 37.0)

# Note that TIME_PERIOD refers to one of our global variables, defined at the top of the script
results = earthaccess.search_data(
    short_name = 'MCD12Q1',
    temporal = TIME_PERIOD,
    bounding_box = tuple(bbox))

# Only download the files once; i.e., if we haven't already downloaded any
if os.path.exists(snakemake.input[0]):
    earthaccess.download(results, os.path.dirname(snakemake.input[0]))

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
