'''
For ET and PET data from VIIRS VNP16 product, clips and projects the
data onto a 1-km global EASE-Grid 2.0.
'''

import datetime
import glob
import json
import os
import numpy as np
import h5py
import rasterio
import rioxarray
import xarray as xr
import py4eos
import pyproj
from rasterio.warp import calculate_default_transform, reproject, Resampling
from tqdm import tqdm

def main():
    # Get a list of all the files
    file_list = glob.glob(f'{snakemake.config["inputs"]["vnp16a2"]}/*.h5')
    file_list.sort()

    et_data = []
    pet_data = []
    for filename in tqdm(file_list):
        # Read the date from the filename
        date = datetime.datetime.strptime(filename.split('/')[-1].split('.')[1], 'A%Y%j')
        # Convert the date to a YYYYMMDD string
        date_str = date.strftime('%Y%m%d')
        # Read in the input VNP16 granule
        hdf = py4eos.read_file(filename, platform = 'VIIRS')
        # Prepare the output filename, project the raster data
        et_filename = f'{snakemake.config["outputs"]["et"].format(date = date_str)}'
        et0 = reproject_viirs(hdf, 'ET_500m', et_filename, driver = 'GTiff')
        pet_filename = f'{snakemake.config["outputs"]["pet"].format(date = date_str)}'
        pet0 = reproject_viirs(hdf, 'PET_500m', pet_filename, driver = 'GTiff')

        # Open the GeoTIFF files as xarray Datasets
        et = rioxarray.open_rasterio(et_filename)\
            .rename(band = 'time')\
            .assign_coords(time = [date])
        et_data.append(et)
        pet = rioxarray.open_rasterio(pet_filename)\
            .rename(band = 'time')\
            .assign_coords(time = [date])
        pet_data.append(pet)

    # Convert from [mm 8day-1] to [mm day-1]
    ds_et0 = xr.concat(et_data, dim = 'time').to_dataset(name = 'ET') / 8.0
    ds_pet0 = xr.concat(pet_data, dim = 'time').to_dataset(name = 'PET') / 8.0
    ds_et = xr.merge([ds_et0, ds_pet0])
    ds_et.to_netcdf(snakemake.output[1])


def reproject_viirs(hdf, field, output_path = '', driver = 'MEM'):
    '''
    Reprojects a VIIRS ET dataset to the global EASE-Grid 2.0.

    Parameters
    ----------
    hdf : py4eos.EOSHDF4
        The EOSHDF4 instance connected to the VIIRS ET dataset
    field : str
        The name of the data variable, e.g., "ET_500m"
    output_path : str
        (Optional) The file path, if writing to a file on disk
    driver : str
        (Optional) The driver name, defaults to "MEM"

    Returns
    -------
    rasterio.io.DatasetWriter
    '''
    et_raster = hdf.to_rasterio(
        field, filename = '', driver = 'MEM', nodata = 32766., scale_and_offset = True)

    # First, resample the ET data to 1-km resolution
    arr = et_raster.read(out_shape = (1200, 1200), resampling = Resampling.average)
    arr = np.where(np.abs(arr) >= 32700, np.nan, arr)
    # We have to re-create the raster dataset, now at 1-km resolution
    et_raster_1km = rasterio.open(
        '', 'w+', driver = 'MEM', height = 1200, width = 1200,
        count = 1, dtype = np.float32, crs = et_raster.crs,
        transform = et_raster.transform * et_raster.transform.scale(2)) # NOTE: Scaling to 1 km
    et_raster_1km.write(arr[0], 1)

    # Second, project the data onto a global EASE-Grid 2.0
    new_transform, width, height = calculate_default_transform(
        et_raster_1km.crs, pyproj.CRS(6933), 1200, 1200, *et_raster_1km.bounds)
    et_raster_ease2 = rasterio.open(
        output_path, 'w+', driver = driver, height = height, width = width,
        count = 1, dtype = np.float32, crs = pyproj.CRS(6933), transform = new_transform)
    reproject(
        source = rasterio.band(et_raster_1km, 1),
        destination = rasterio.band(et_raster_ease2, 1),
        resampling = Resampling.bilinear,
        src_nodata = np.nan, # Necessary so that missing data is interpolated
        dst_nodata = np.nan)
    return et_raster_ease2



if __name__ == '__main__':
    main()
