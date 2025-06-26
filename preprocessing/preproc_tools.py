# IMPORT LIBRARIES

import xarray as xr
#import rioxarray as rio
from rasterio.warp import Resampling
import numpy as np
#from scipy.ndimage import convolve
import pandas as pd
import geopandas as gpd
import os
from shapely.geometry import mapping
import glob
import time
import requests
#from itertools import product
from zipfile import ZipFile

## Download functions

def download_spam(refyear, spam_variable, basepath):
    # Prepare download url and path
    if refyear == '2010':
        if spam_variable == 'physical_area':
            url = "https://s3.amazonaws.com/mapspam-data/2010/v2.0/geotiff/spam2010v2r0_global_phys_area.geotiff.zip"
        if spam_variable == 'harvested_area':
            url = "https://s3.amazonaws.com/mapspam-data/2010/v2.0/geotiff/spam2010v2r0_global_harv_area.geotiff.zip"
        if spam_variable == 'yield':
            url = "https://s3.amazonaws.com/mapspam/2010/v2.0/geotiff/spam2010v2r0_global_yield.geotiff.zip"
        if spam_variable == 'production':
            url = "https://s3.amazonaws.com/mapspam-data/2010/v2.0/geotiff/spam2010v2r0_global_prod.geotiff.zip"
    if refyear == '2020': # DOES NOT WORK YET
        if spam_variable == 'physical_area':
            url = "https://www.dropbox.com/scl/fi/napqtql4521ujqt22j05w/spam2020V1r0_global_physical_area.geotiff.zip?rlkey=vpamm4zj3gu2752ubpj3j80iu&dl=0"
        if spam_variable == 'yield':
            url = ""    # TO BE UPDATED WHEN SPAM FIXES THE WEBLINK

    if (spam_variable == 'physical_area') or (spam_variable == 'harvested_area'):
        target_dir = basepath + 'preprocessing/rawdata/cropmasks/'
    if (spam_variable == 'yield') or (spam_variable == 'production'):
        target_dir = basepath + 'preprocessing/rawdata/calibration/'
    if not os.path.exists(target_dir):
        os.mkdir(target_dir)
    download_path = target_dir + 'spam' + refyear + '_' + spam_variable + '.zip'

    # Start download
    print("        *** DOWNLOADING SPAM " + refyear + " data (" + spam_variable + ") ***")
    print('        URL: ', url)
    download_url(url, download_path=download_path)
    print("        *** UNZIPPING SPAM " + refyear + " data (" + spam_variable + ") ***")
    unzip_all(dir=target_dir)
    unzipped_download_directory = download_path[:-4]

    # Remove unnecessary files -> Rather delete entire raw data folder at the end of preprocessing!
    #for filename in glob.glob(download_path[:-3]+"target_dir/version*"):
    #    os.remove(filename)

    return unzipped_download_directory

def download_url(url, download_path="./"):
    """Download data from a given url to a given download path on the hard drive"""
    while True:
        try:
            response = requests.get(url, timeout=30)
        except:
            print("    *** WARNING! Download froze, likely due to network instability. Restarting download...")
            time.sleep(2)
        else:
            if response.ok:
                # inmemory = tiff.imread(io.BytesIO(response.content))  # reads download data directly to memory (without writing to disk)
                with open(download_path, 'wb') as savefile:
                    savefile.write(response.content)  # write download data to disk
                break

def unzip_all(dir='.'):
    """Iteratively unzip files in given directory and subdirectories until no zip files remain"""
    zipfiles = True
    while bool(zipfiles):
        zipfiles = []
        for root, d_names, f_names in os.walk(dir):
            for f in f_names:
                if f.endswith('.zip'):
                    zipfiles.append(os.path.join(root, f))
        if bool(zipfiles):
            for file in zipfiles:
                # targetdir = dir + '/' + file[:-4]
                targetdir = file[:-4]
                if not os.path.exists(targetdir):
                    os.mkdir(targetdir)
                with ZipFile(file, 'r') as zObject:
                    zObject.extractall(path=targetdir)
                os.remove(file)

## Preprocessing functions
def preproc_spam(basepath, download_dir, refyear, spam_variable, domain_path, referenceraster_path):
    mask = gpd.read_file(domain_path)
    to_match = xr.load_dataset(referenceraster_path)
    to_match.rio.write_crs(4326, inplace=True)

    # Dictionary that connects FAO crop ID's to crop names in AquaCrop. Used for layer naming
    crop_dict = {'BARL': 'Barley','COTT': 'Cotton','BEAN': 'DryBean','MAIZ': 'Maize','RICE': 'PaddyRice','POTA': 'Potato','SORG': 'Sorghum','SOYB': 'Soybean','SUGB': 'SugarBeet','SUGC': 'SugarCane','SUNF': 'Sunflower','WHEA': 'Wheat_summer','CASS': 'Cassava'}

    # Dictionnary to rename technology
    tech_dict = {'R': 'rf', 'I': 'ir'}

    # CHECK IF REALLY ALL FOUR VARIABLES ARE NEEDED. PROBABLY JUST PHYSICAL AREA AND - IF CALBRATION INCLUDED - YIELD
    search_criteria = download_dir + '/*.tif'
    layers = glob.glob(search_criteria)
    file_to_mosaic = []
    for layer in layers:

        # Preserve only layers referring to rainfed (R) or irrigated (I) technology
        technique = layer[-5:-4]
        if technique not in ['R','I']:
            os.remove(layer)
            continue

        # Get crop ID from layer name and preserve only if crop type is supported by AquaCrop
        cropID = layer[-10:-6]
        if cropID not in crop_dict:
            os.remove(layer)
            continue

        print(cropID, technique)

        src = xr.load_dataset(layer)
        src_proj = src.rio.reproject_match(to_match, resampling=Resampling.average) # changed from nearest-neighbor to area-weighted resampling to reduce distortions
        src_clip = src_proj.rio.clip(mask.geometry.apply(mapping))

        # Rename variables from FAO crop ID to crop strings used by AquaCrop
        layername_base = crop_dict.get(cropID) + '_' + tech_dict.get(technique)

        # Change resolution from 5 arcmin to 3 arcmin areas and percentages for physical and harvested area and production rasters
        if (spam_variable == 'physical_area') | (spam_variable == 'harvested_area') | (spam_variable == 'production'):
            src_clip[layername_base + '_' + spam_variable] = src_clip['band_data'] * 3**2 / 5**2
            #src_clip[layer[-12: -4] + '_percentage'] = src_clip[layer[-12: -4] + '_area'] / meters.ha  # Percentage only needed to restrict model to cells with very small crop area

        # Renaming variable for yield and change from kg/ha to t/ha
        elif spam_variable == 'yield':
            src_clip[layername_base + '_yield'] = src_clip['band_data'] / 1000

        # Dropping band_data for merging
        src_clip = src_clip.drop_vars(['band_data'])

        file_to_mosaic.append(src_clip)

    # Merge data into one file
    src_mosaic = xr.merge(file_to_mosaic)
    src_mosaic.to_netcdf(basepath + 'preprocessing/processed/' + 'spam' + refyear + '_' + spam_variable + '.nc')

def spam_refyear(start_year, end_year):
    # Choose SPAM data reference year (available reference years are 2000, 2005, 2010, or 2020) according to average of start and end year of modelling horizon
    import numpy as np
    avg_year = np.mean([start_year, end_year])
    avg_year = np.ceil(avg_year)
    SPAM_refyears = [2000, 2005, 2010, 2020]
    refyear = min(SPAM_refyears, key=lambda x: abs(x - avg_year)) # select nearest reference year available in SPAM

    return refyear

## Preprocessing climate data
def preproc_era5(src, variable, yearlist, basepath, to_match):
    print("        *** PREPROCESSING CLIMATE DATA: " + variable + " ***")

    # Preparations
    src = ensure_xy_dims(src)
    src.rio.write_crs(4326, inplace=True)

    varname_dict = {'MinTemp': 't2m', 'MaxTemp': 't2m', 'Precipitation': 'tp', 'ReferenceET': 'pev'}  # Names of data variables in AquaCrop and ERA5 data, respectively

    # Adjust units and data variable names to AquaCrop definitions
    if variable in ['MinTemp', 'MaxTemp']:
        src = src.rename_vars({varname_dict.get(variable): variable})  # Rename data variables to AquaCrop names
        src = src - 273.15  # Convert from Kelvin to Celsius
    if variable in ['Precipitation', 'ReferenceET']:
        src = src.rename_vars({varname_dict.get(variable): variable})  # Rename data variables to AquaCrop names
        src = src * 1000 # Convert from metres to millimetres
        if variable in ['ReferenceET']:
            src = -src  # Multiply by -1 since raw data is negative due to ERA5 definition (i.e. downward fluxes are positive, upward fluxes negative)
        elif variable in ['Precipitation']:
            # Adjust time coordinate to date of previous day, as Jan 1 00:00 represents daily accumulation of Dec 31 of previous year (see ERA5 time step definition: https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation )
            src = src.isel(valid_time=slice(1, None))   # Delete first time step as it represents daily accumulation of 31/12 of previous year
            valid_dates = pd.to_datetime(src.valid_time.values)  # Convert valid_time to pandas DatetimeIndex
            new_dates = (valid_dates.normalize() - pd.Timedelta(days=1))  # Strip time part (keep just the date) and shift back by one day
            src = src.assign_coords(valid_time=new_dates)  # Assign adjusted timestamps back to the dataset
            src = src.drop_vars('expver')
        src[variable] = src[variable].where(src[variable] >= 0, 0)  # Set negative precipitation and evaporation values to 0.

    # Fill missing values in data variable(s)
    import scipy.ndimage as ndi
    variables = list(src.data_vars)
    for variab in variables:    # Loop to make it work also for InitSoilwater data, which has four data variables (corresponding to four soil depths)
        data3d = src[variab].to_numpy()
        nodata_mask = np.isnan(data3d[0])  # This assumes that all timesteps have the same missing cells (fair assumption, since the same ERA5-land water mask is applied to all time steps)
        dist, nearest_indices = ndi.distance_transform_edt(nodata_mask, return_indices=True)  # Nearest-neighbour interpolation
        src[variab].data = data3d[:, nearest_indices[0], nearest_indices[1]]

    # Resample to project grid and mask
    src_reproj = src.rio.reproject_match(to_match, resampling=Resampling.nearest)
    src_reproj = src_reproj.rename({'valid_time': 'time', 'x': 'longitude', 'y': 'latitude'})
    src_masked = src_reproj.where(to_match['Band1'] == 1)

    #if variable in ['InitSoilwater']:  # NOT NEEDED AFTER ALL BECAUSE AQUACROP CAN DO IT USING ORIGINAL SOIL LAYERS FROM ERA5 (see function get_initial_WC in aquacropgrid.py)
    #    src_reproj = convert_soildepthlayers(src_reproj)

    # Prepare output directory
    target_dir = basepath + 'preprocessing/processed/'
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    # Save to disk
    src_masked.to_netcdf(target_dir + variable + str(yearlist[0]) + str(yearlist[-1]) + '.nc', mode='w', encoding = {variables[0]: {'zlib': True, 'complevel': 4}})  # Save to disk with moderate compression

def convert_soildepthlayers(src_in):    # PROBABLY REMOVE AS NOT NEEDED
    """
    Take soil water content data from four soil depth layers (as defined in ERA5-Land) and convert it into the six soil layers defined by AquaCrop.
    We simply match the layers here as follows:
    ERA5-Land layer 1: 0-7cm            AquaCrop layer 1:         0-5cm
    ERA5-Land layer 2: 7-28cm           AquaCrop layers 2 & 3:    5-15cm, 15-30cm
    ERA5-Land layer 3: 28-100cm         AquaCrop layers 4 & 5:    30-60cm, 60-100cm
    ERA5-Land layer 4: 100-289cm        AquaCrop layer 6:         100-200cm
    """
    src = src_in.copy()

    # Rename data variables according to AquaCrop soil layer definition
    rename_dict = {'swvl1': 'initwc_layer1', 'swvl2': 'initwc_layer2', 'swvl3': 'initwc_layer4', 'swvl4': 'initwc_layer6'}
    src = src.rename_vars(rename_dict)

    # Add copies of data for AquaCrop soil layers 3 and 5
    src['initwc_layer3'] = src['initwc_layer2']
    src['initwc_layer5'] = src['initwc_layer4']

    # Reorder the dataset
    order = ['initwc_layer1', 'initwc_layer2', 'initwc_layer3', 'initwc_layer4', 'initwc_layer5', 'initwc_layer6']
    src_ordered = src[order]

    return src_ordered


## Helper functions

def basegrid(domain_shape_path, resolution, templategrid_path):    # Creates basic raster file as template for all preprocessing scripts (i.e. what is read in all other scripts as "to_match" file)
    from rasterio.transform import from_origin
    from rasterio import features
    from affine import Affine

    # Read shapefile
    mask = gpd.read_file(domain_shape_path)

    # Check if all geometries are Polygon or MultiPolygon
    if not all(g in ['Polygon', 'MultiPolygon'] for g in mask.geom_type.unique()):
        raise Exception("The input polygon file contains geometries other than Polygon/MultiPolygon.")

    # Check if it's EPSG:4326
    if not mask.crs.to_epsg() == 4326:
        raise Exception("The polygon file must be in EPSG:4326.")

    full_geom = mask.unary_union
    xmin, ymin, xmax, ymax = full_geom.bounds # in lat/lon
    xmin, ymin = np.floor(xmin * (1/resolution)) / (1/resolution) , np.floor(ymin * (1/resolution)) / (1/resolution) # round bounds down to cell resolution
    xmax, ymax = np.ceil(xmax * (1/resolution)) / (1/resolution) , np.ceil(ymax * (1/resolution)) / (1/resolution) # round upper bounds up to cell resolution
    bounds = [xmin, ymin, xmax, ymax]

    if os.path.exists(templategrid_path):
        ds = xr.open_dataset(templategrid_path)

    else:
        w = round((xmax - xmin) / resolution)
        h = round((ymax - ymin) / resolution)
        lons = np.arange(xmin+resolution/2, xmax, resolution)
        lats = np.arange(ymin+resolution/2, ymax, resolution)

        # Affine transform: top-left origin
        transform = Affine.translation(xmin, ymax) * Affine.scale(resolution, -resolution)

        # Mask: True = inside polygon
        mask = features.geometry_mask([full_geom],out_shape=(h, w),transform=transform,invert=True)

        # Create data array: inside = 1, outside = nodata
        nodata_value = 0
        data = np.full((h, w), nodata_value, dtype=np.uint8)
        data[mask] = 1

        # CF convention prefers descending latitudes
        lats = lats[::-1]

        # Create xarray Dataset
        ds = xr.Dataset({"Band1": (["latitude", "longitude"], data),},coords={"latitude": lats,"longitude": lons,},attrs={"Conventions": "CF-1.8","title": "Template raster from polygon shapefile","crs": "EPSG:4326"})
        ds["spatial_ref"] = xr.DataArray(0, attrs={"grid_mapping_name": "latitude_longitude", "epsg_code": 4326, "semi_major_axis": 6378137.0, "inverse_flattening": 298.257223563, "long_name": "CRS definition"})
        ds["Band1"].attrs.update({"grid_mapping": "spatial_ref", "_FillValue": nodata_value, "missing_value": nodata_value})

        # Save to NetCDF
        ds.to_netcdf(templategrid_path)

    return ds, bounds

# Ensures the dataset has 'x' and 'y' dimensions (for resampling in xarray)
def ensure_xy_dims(ds):
    # Rename dimensions to 'x' and 'y' if they are named differently
    if 'longitude' in ds.dims and 'latitude' in ds.dims:
        ds = ds.rename({'longitude': 'x', 'latitude': 'y'})
    elif 'lon' in ds.dims and 'lat' in ds.dims:
        ds = ds.rename({'lon': 'x', 'lat': 'y'})
    return ds

