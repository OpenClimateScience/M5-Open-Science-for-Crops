'''
Downloads ET and PET data from VIIRS VNP16 product, then clips and projects the
data onto a 1-km global EASE-Grid 2.0.
'''

import datetime
import glob
import os
import numpy as np
import h5py
import rasterio
import py4eos
import pyproj
from rasterio.warp import calculate_default_transform, reproject, Resampling
from tqdm import tqdm

def main():
    # Get a list of all the files
    file_list = glob.glob(f'{snakemake.config["paths"]["dir_vnp16a2"]}/*.h5')
    file_list.sort()

    for filename in tqdm(file_list):
        # Read the date from the filename
        date = datetime.datetime.strptime(filename.split('/')[-1].split('.')[1], 'A%Y%j')
        # Convert the date to a YYYYMMDD string
        date_str = date.strftime('%Y%m%d')
        # Prepare the output filename
        output_file_tpl = f'{snakemake.config["paths"]["dir_vnp16a2"]}/VNP16_%s_mm_8day-1_{date_str}.tiff'
        hdf = py4eos.read_file(filename, platform = 'VIIRS')
        et = reproject_viirs(hdf, 'ET_500m', output_file_tpl % 'ET', driver = 'GTiff')
        pet = reproject_viirs(hdf, 'PET_500m', output_file_tpl % 'PET', driver = 'GTiff')


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
