'''
Processes the IMERG-Final precipitation data into a data cube, resampled onto
a 1-km global EASE-Grid 2.0.
'''

import datetime
import glob
import os
import numpy as np
import h5py
import rasterio
import py4eos
import pyproj
from shapely.geometry import Polygon
from rasterio.warp import calculate_default_transform, reproject, Resampling
from tqdm import tqdm

et_raster_ease2 = rasterio.open(snakemake.input[0])

bb = et_raster_ease2.bounds
bounds = Polygon([
    (bb.left, bb.bottom),
    (bb.left, bb.top),
    (bb.right, bb.top),
    (bb.right, bb.bottom)
])

file_list = glob.glob(f"{snakemake.config['paths']['dir_imerg']}/*.nc4")
file_list.sort()

stack = []
for filename in tqdm(file_list):
    ds = xr.open_dataset(filename)
    # Anytime we use the .rio attribute, we must have rioxarray installed
    ds_ease2 = ds[['precipitation']]\
        .transpose('time', 'lat', 'lon')\
        .rio.write_crs(4326)\
        .rio.set_spatial_dims('lon', 'lat')\
        .rio.reproject(pyproj.CRS(6933), resolution = 9000)\
        .rio.clip([bounds])
    stack.append(ds_ease2)
