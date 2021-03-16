# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 14:41:40 2019
Calculate direct flux gain and loss by ecoregion
@author: pj276
"""

#%%
# Imports
import geopandas as gpd, rasterio, os, ogr, sys
from shapely.geometry import box
import fiona, gdal
import pandas as pd
from rasterio.tools.mask import mask
import numpy as np
sys.path.append('/home/pj276/projects/ifl_corridors/code/functions/')
import hf_utilities as hfu
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
import param_file as p1

#%%
# Top level dir
idir = p1.idir #'/scratch/pj276/ifl_corridors'

# Directory that will hold output direct flux files
fluxdir = p1.fluxdir

#%%
# Get rastproj
rtemp = p1.rtemp #'/scratch/pj276/ifl_corridors/rp5k130k90m_neg01_as_2013_ifl/cfinal/flux/corr_2013_120_99_flux.tif'
rtemp = rasterio.open(rtemp)
rproj4 = rtemp.crs

#%%
# Read in zone shapefile
inshp = p1.inshp 
shpfile = gpd.read_file(inshp, encoding="utf-8")

# Specify field
field = p1.field

# Get list of unique names/codes
ucode = np.unique(shpfile[field].tolist())

# Job array index
jai = int(sys.argv[1]) - 1 # subtract for python zero indexing

#%%
# Region
i = p1.region #['af','as','sa']
# Scenario
j = p1.scenario #['neg2', 'neg01', 'pos2']
# Year
year = ['2000', '2013', '2018']

# Empty list for shapefile names
sflist = []

# Loop through and calculate flux by zone
for k in year:
    # We need kk to allow the directory year to be different from the file year
    # since the 2018 loss refers to the 2013 flux layer
    # Name of direct flux raster
    if k == '2000':
        kk = '2000'
        dfrname = idir + '/rp5k130k90m_' + j + '_' + i + '_' + kk + '_ifl/cfinal/flux/dfloss_rp5k130k90m_' + j + '_' + i + '_' + kk + '_mos.tif'
    if k == '2013':
        kk = '2013'
        dfrname = idir + '/rp5k130k90m_' + j + '_' + i + '_' + kk + '_ifl/cfinal/flux/dfgain_rp5k130k90m_' + j + '_' + i + '_' + kk + '_mos.tif'
    if k == '2018':
        kk = '2013'
        dfrname = idir + '/rp5k130k90m_' + j + '_' + i + '_' + kk + '_ifl/cfinal/flux/dfloss1318_rp5k130k90m_' + j + '_' + i + '_' + kk + '_mos.tif'
    # Get direct flux raster corners
    minx,miny,maxx,maxy = hfu.raster_corners(dfrname)
    b = [box(minx, miny, maxx, maxy)]
    d = {'id': [1], 'region': [i]}
    dd = pd.DataFrame(data=d)
    bbox = gpd.GeoDataFrame(dd, geometry=b)
    bbox.crs = rproj4
    bboxgeog = bbox.to_crs(shpfile.crs) 

    # Test if unique code ucode is string or unicode, treat differently if so
    if isinstance(ucode[jai], str):
        l = ucode[jai]
        lorig = l
        if ' ' in l:
            l = l.replace(' ', '_')
    elif isinstance(ucode[jai], unicode):
        l = ucode[jai]
        lorig = l
        if ' ' in l:
            l = l.replace(' ', '_')
    else:
        l = ucode[jai]
        lorig = l
        l = np.str(np.int(l))
    # Output tiff name
    otiff1 = idir + '/rp5k130k90m_' + j + '_' + i + '_' + kk + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '_dflux.tif'

    # Get shapes
    cx = shpfile.loc[[s for s, n in enumerate(shpfile[field]) if n == lorig],]
    # If at least one shape is returned, keep processing
    if cx.shape[0] > 0:
        # Get target column and geometry column
        cx = cx[[field,'geometry']]
        # Dissolve if there are multiple shapes
        if cx.shape[0] > 1:
            cx = cx.dissolve(by=field, aggfunc='first')
            cx[field] = cx.index
        # Fix simple self intersections if necessary
        if cx.is_valid.bool() == False:
            cx = cx.buffer(0)
        # Check if polygon intersects with raster
        minx, miny, maxx, maxy = cx.geometry.total_bounds
        b = [box(minx, miny, maxx, maxy)]
        d = {'id': [1], 'ucode': [l]}
        dd = pd.DataFrame(data=d)
        bboxcx = gpd.GeoDataFrame(dd, geometry=b)
        bboxcx.crs = cx.crs
        ix = gpd.overlay(bboxgeog, bboxcx, how='intersection', use_sindex=True)
        if ix.shape[0] >= 1:
            # Repro shapefile
            cx = cx.to_crs(rproj4)
            # Save shapefile if needed
            if not os.path.exists(idir + '/rp5k130k90m_' + j + '_' + i + '_' + kk + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '.shp'):
                cx.to_file(idir + '/rp5k130k90m_' + j + '_' + i + '_' + kk + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '.shp')              
                sflist.append(idir + '/rp5k130k90m_' + j + '_' + i + '_' + kk + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '.shp')
            # Get shapefile extent in regional sinusoidal coordinates
            minx, miny, maxx, maxy = cx.geometry.total_bounds
            # Read in direct flux raster, subset by shape bounding box
            dfr = hfu.raster2array(dfrname, clip=[minx, miny, maxx, maxy])
            # If raster2array returns a valid output, keep processing
            if dfr != False:
                # Write to file
                rwkt = gdal.Open(dfrname)
                rwkt = rwkt.GetProjectionRef()
                hfu.array2raster(otiff1,'#',gdal.GDT_Float32,'GTiff',dfr[0],geotrans=dfr[1],rasterproj=rwkt)
                # Get geometry from zone shapefile
                with fiona.open(idir + '/rp5k130k90m_' + j + '_' + i + '_' + kk + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '.shp', "r") as shapefile:
                    geoms = [feature["geometry"] for feature in shapefile]
                # Use zone boundary to mask flux sum layer
                with rasterio.open(otiff1) as src:
                    out_image, out_transform = mask(src, geoms, crop=True)
                    out_meta = src.meta.copy()
                # Delete rectangle cropped raster
                os.remove(otiff1)
                # Sum flux values to get total for the zone
                tcsum = np.sum(out_image)
                # One line dataframe to hold flux and other info
                d = {'region': [i], 'scenario': [j], 'year': [k], 'ucode': [l], 'flux': [tcsum]}
                fluxdf = pd.DataFrame(data=d)
                # Write to file
                fluxdf.to_csv(idir + '/rp5k130k90m_' + j + '_' + i + '_' + kk + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '_flux.csv')
# Delete shapefiles
driver = ogr.GetDriverByName("ESRI Shapefile")
for p in list(set(sflist)):
    driver.DeleteDataSource(p)
