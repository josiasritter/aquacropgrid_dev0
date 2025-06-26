## Import packages and prepare

import xarray as xr
#import rioxarray as rio
import numpy as np
#from rasterio.warp import Resampling
import os
import geopandas as gpd
import cdsapi
from shapely.geometry import mapping
import sys
sys.path.append("./preprocessing")
from preproc_tools import preproc_era5, basegrid

"""
Download from ERA5 / ERA5-Land and prepocess: 
    - Four daily climate variables: 
        - Minimum temperature [K]
        - Maximum temperature [K]
        - Precipitation [m]
        - Potential evapotranspiration [m]
    - Soil moisture of initial time step [m3/m3]
"""

## INPUT ARGUMENTS. REPLACE THESE WITH YOUR OWN CONFIGURATION (for testing; will be replaced by function inputs)
#def climate(basepath, domain_path, referenceraster_path, start_year, end_year, api_token, cell_resolution, variables=['temperature_daily_minimum','temperature_daily_maximum','precipitation','potevaporation','InitSoilwater']):

start_year = 2020
end_year = 2021
api_token = 'XXX'  # your API token, retrieved from your profile page on the Copernicus Climate Data Store (https://cds.climate.copernicus.eu/)
basepath = os.path.abspath(os.path.join(os.path.dirname(__file__), '..')) + '/'   # your home directory
domain_path = os.path.join(basepath, 'inputdata', 'mekong', 'basin_outline', 'mekong_jrc_outline.geojson')     # this polygon file defines the area of interest (must be in lat/lon, i.e. EPSG 4326)
cell_resolution = 0.05    # degrees
variables=['MinTemp','MaxTemp','Precipitation','ReferenceET','InitSoilwater']    # testing all variables
#variables=['Precipitation'] # testing individual variables

## Years to be downloaded
yearlist = list(range(start_year, end_year+1))

# Area to be downloaded (bounding box)
templategrid_path = basepath + 'preprocessing/template_grid.nc'
to_match, bounds = basegrid(domain_path, cell_resolution, templategrid_path)
#to_match = xr.open_dataset(templategrid_path)
#to_match.rio.write_crs(4326, inplace=True)
#boundbox = to_match.rio.bounds()
#bounds = np.around(bounds, decimals=2).tolist()
bounds=[bounds[3],bounds[0],bounds[1],bounds[2]]    # reorder bounds to follow ERA5 CDS definition (N-W-S-E)

# Prepare download directory
target_dir = basepath + 'preprocessing/rawdata/climate/'
if not os.path.exists(target_dir):
    os.makedirs(target_dir)

# Prepare Copernicus Climate Data Store (CDS) API
#os.environ['CDSAPI_URL'] = 'https://cds.climate.copernicus.eu/api'
url = 'https://cds.climate.copernicus.eu/api'
c = cdsapi.Client(url=url, key=api_token)

## Download daily min and max temperatures from ERA5-Land daily

t_stats = ["daily_minimum", "daily_maximum"]
for stat in t_stats:
    variable = "MinTemp" if stat == "daily_minimum" else "MaxTemp"  # rename variables to AquaCrop definition
    if variable in variables:
        for year in yearlist:   # Split download into separate years
            if not os.path.exists(target_dir + variable + str(year) + '.nc'):   # Skip download if file already exists
                print("        *** DOWNLOADING CLIMATE DATA: " + variable + str(year) + " ***")
                c.retrieve(
                    'derived-era5-land-daily-statistics',
                    {
                        'variable': ['2m_temperature'],
                        'year': [str(year)],
                        'month': ["01","02","03","04","05","06","07","08","09","10","11","12"],
                        "day": ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"],
                        "daily_statistic": stat,
                        "time_zone": "utc+00:00",
                        "frequency": "6_hourly",
                        "area": bounds,
                        'data_format': 'netcdf',
                    },
                    target_dir + variable + str(year) + '.nc'
                )

        # Combine yearly files into one
        file_paths = [target_dir + variable + str(year) + '.nc' for year in yearlist]
        #src = xr.open_mfdataset(file_paths, combine='nested', concat_dim='valid_time', engine='netcdf4', chunks={})    # Fails because original files are GRIB files that are not compatible with dask processing, even when converted to netcdf
        datasets = [xr.open_dataset(f) for f in file_paths]
        src = xr.concat(datasets, dim='valid_time')
        src = src.sortby('valid_time')    # Make sure all is properly sorted along the time dimension

        # Preprocessing
        preproc_era5(src, variable, yearlist, basepath, to_match)

## Download daily precipitation accumulation from ERA5-Land hourly (time step 0000 UTC represents full accumulation of previous day, see https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation)

variable = 'Precipitation'
if variable in variables:
    for year in yearlist:
        if not os.path.exists(target_dir + variable + str(year) + '.nc'):   # Skip download if file already exists
            print("        *** DOWNLOADING CLIMATE DATA: " + variable + str(year) + " ***")
            c.retrieve(
                "reanalysis-era5-land",
                {
                    "variable": ["total_precipitation"],
                    "year": [str(year)],
                    "month": ["01","02","03","04","05","06","07","08","09","10","11","12"],
                    "day": ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"],
                    "time": ["00:00"],
                    "area": bounds,
                    "data_format": "netcdf",
                    "download_format": "unarchived"
                },
                target_dir + variable + str(year) + '.nc'
            )

    # Last time step (Dec 31) is saved in Jan 1 00:00 of following year -> Separate download necessary with this (annoying) ERA5 definition of accumulation variables in hourly data
    if not os.path.exists(target_dir + variable + str(end_year+1) + 'Jan1.nc'):  # Skip download if file already exists
        print("        *** DOWNLOADING CLIMATE DATA: " + variable + str(end_year+1) + "Jan1 ***")
        c.retrieve(
            "reanalysis-era5-land",
            {
                "variable": ["total_precipitation"],
                "year": [str(end_year+1)],
                "month": ["01"],
                "day": ["01"],
                "time": ["00:00"],
                "area": bounds,
                "data_format": "netcdf",
                "download_format": "unarchived"
            },
            target_dir + variable + str(end_year+1) + 'Jan1.nc'
        )

    # Combine yearly files into one
    file_paths = [target_dir + variable + str(year) + '.nc' for year in yearlist] + [target_dir + variable + str(end_year+1) + 'Jan1.nc']
    datasets = [xr.open_dataset(f) for f in file_paths]
    src = xr.concat(datasets, dim='valid_time')
    src = src.sortby('valid_time')  # Make sure all is properly sorted along the time dimension

    # Preprocessing
    preproc_era5(src, variable, yearlist, basepath, to_match)

## Download daily potential evaporation accumulation from ERA5 daily (ERA5 used since ERA5-Land defines potential evaporation differently)
variable = 'ReferenceET'
if variable in variables:
    for year in yearlist:
        if not os.path.exists(target_dir + variable + str(year) + '.nc'):  # Skip download if file already exists
            print("        *** DOWNLOADING CLIMATE DATA: " + variable + str(year) + " ***")
            c.retrieve(
                "derived-era5-single-levels-daily-statistics",
                {
                    "product_type": "reanalysis",
                    "variable": ["potential_evaporation"],
                    "year": [str(year)],
                    "month": ["01","02","03","04","05","06","07","08","09","10","11","12"],
                    "day": ["01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"],
                    "daily_statistic": "daily_sum",
                    "time_zone": "utc+00:00",
                    "frequency": "1_hourly",
                    "area": bounds,
                    "data_format": "netcdf"
                },
                target_dir + variable + str(year) + '.nc'
            )

    # Combine yearly files into one
    file_paths = [target_dir + variable + str(year) + '.nc' for year in yearlist]
    datasets = [xr.open_dataset(f) for f in file_paths]
    src = xr.concat(datasets, dim='valid_time')
    src = src.sortby('valid_time')    # Make sure all is properly sorted along the time dimension

    # Preprocessing
    preproc_era5(src, variable, yearlist, basepath, to_match)

## Download Volumetric soil water content [m3/m3] for initial time step (has four soil depth layers) from ERA5-Land hourly

variable = 'InitSoilwater'
if variable in variables:
    if not os.path.exists(target_dir + variable + str(start_year) + '.nc'):  # Skip download if file already exists
        print("        *** DOWNLOADING CLIMATE DATA: " + variable + " ***")
        c.retrieve(
            "reanalysis-era5-land",
            {
                "variable": [
                    "volumetric_soil_water_layer_1",
                    "volumetric_soil_water_layer_2",
                    "volumetric_soil_water_layer_3",
                    "volumetric_soil_water_layer_4"
                ],
                "year": [str(start_year)],
                "month": ["01"],
                "day": ["01"],
                "time": ["00:00"],
                "data_format": "netcdf",
                "download_format": "unarchived",
                "area": bounds
            },
            target_dir + variable + str(start_year) + '.nc'
        )

    # Preprocessing
    src = xr.open_dataset(target_dir + variable + str(start_year) + '.nc')
    preproc_era5(src, variable, yearlist, basepath, to_match)
