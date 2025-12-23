'''
Downloads the MCD12Q1 land-cover data.
'''

import glob
import os
import earthaccess

auth = earthaccess.login()

TIME_PERIOD = ('2023-10-01', '2024-09-30')

# The bounding box for our study area
bbox = (1.5, 34.0, 8.0, 37.0)

results_mcd12q1 = earthaccess.search_data(
    short_name = 'MCD12Q1',
    temporal = TIME_PERIOD,
    bounding_box = tuple(bbox))
earthaccess.download(results_mcd12q1, os.path.dirname(snakemake.params.dir_mcd12q1))
