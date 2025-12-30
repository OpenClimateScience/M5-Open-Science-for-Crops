'''
Processes the IMERG-Final precipitation data into a data cube, resampled onto
a 1-km global EASE-Grid 2.0.
'''

import glob
import json
import os
import numpy as np
import pyproj
import xarray as xr
import shapely
from rasterio.warp import calculate_default_transform, reproject, Resampling
from tqdm import tqdm

file_list = glob.glob(f"{snakemake.config['inputs']['imerg']}/*.nc4")
file_list.sort()

with open(snakemake.input[0], 'r') as file:
    bounds = shapely.from_geojson(json.load(file))

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

ds_precip = xr.concat(stack, dim = 'time')
ds_precip.to_netcdf(snakemake.output[0])
