# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 11:12:50 2019
Calculate forest cover in corridors by country
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
shpfile = gpd.read_file(inshp)

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
year = ['2000', '2013']
# No data value for forest cover
ndval = -9999

# Loop through and calculate flux by zone
for k in year:

    # Name of forest cover raster (non-corridor areas have a value of -9999)
    fcname = idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/flux/fc_rp5k130k90m_' + j + '_' + i + '_' + k + '_mos.tif'
    # Get forest cover raster corners
    minx,miny,maxx,maxy = hfu.raster_corners(fcname) 
    bbox = gpd.GeoDataFrame({'id': [1], 'region': [i]}, geometry=[box(minx, miny, maxx, maxy)])
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
    otiff1 = idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '_fcincorrnoifl.tif'

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
        bboxcx = gpd.GeoDataFrame({'id': [1], 'ucode': [l]}, geometry=[box(minx, miny, maxx, maxy)])
        bboxcx.crs = cx.crs
        ix = gpd.overlay(bboxgeog, bboxcx, how='intersection', use_sindex=True)
        if ix.shape[0] >= 1:
            # Save file if needed
            cx = cx.to_crs(rproj4)
            cx.to_file(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '.shp')              
            # Get shapefile extent in regional sinusoidal coordinates
            minx, miny, maxx, maxy = cx.geometry.total_bounds
            # Convert to bounding box
            bboxcx = gpd.GeoDataFrame({'id': [1], 'ucode': [l]}, geometry=[box(minx, miny, maxx, maxy)])
            # Get intersection between raster bounding box and shapefile bounding box in sinusoidal
            ix = gpd.overlay(bbox, bboxcx, how='intersection', use_sindex=True)
            # Get clip coordinates from intersection
            minx, miny, maxx, maxy = ix.geometry.total_bounds
            # Read in forest cover raster, subset by intersection between polygon
            # bounding box and raster.
            dfr = hfu.raster2array(fcname, clip=[minx, miny, maxx, maxy])
            # If raster2array returns a valid output, keep processing
            if dfr != False:
                # Rasterized ifl name
                rifl = '/scratch/pj276/ifl_corridors/ifl/ifl_2000_tropics_ssu_' + i + '_' + j + '.tif'
                # Read in rasterized ifl, subset by shape bounding box
                rifl = hfu.raster2array(rifl, clip=[minx, miny, maxx, maxy])[0]
                # Mask out forest cover that overlaps ifls
                dfr[0][rifl >= 1] = ndval
                dfra = dfr[0]
                # Write to file
                rwkt = gdal.Open(fcname)
                rwkt = rwkt.GetProjectionRef()
                hfu.array2raster(otiff1,'#',gdal.GDT_Float32,'GTiff',dfra,geotrans=dfr[1],rasterproj=rwkt)
                # Get geometry from zone shapefile
                with fiona.open(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '.shp', "r") as shapefile:
                    geoms = [feature["geometry"] for feature in shapefile]
                # Use zone boundary to mask forest cover  layer
                with rasterio.open(otiff1) as src:
                    out_image, out_transform = mask(src, geoms, crop=True, nodata=ndval)
                    out_meta = src.meta.copy()
                # Delete rectangle cropped raster
                os.remove(otiff1)
                # Forest cover area
                fcasum = np.sum(out_image[out_image != ndval]/100.0*8100.0)
                # Stop if there's no tree cover in the zone
                if (fcasum <= 0) or not (np.isfinite(fcasum)):
                    print('invalid or zero tree cover in zone, exiting')
                    sys.exit(0)
                # Average forest cover values for zone
                fcave = np.mean(out_image[out_image != ndval])
                # Standard deviation of forest cover values for zone
                fcstd = np.std(out_image[out_image != ndval])
                # One line dataframe to hold forest cover mean, std, and other info
                d = {'region': [i], 'scenario': [j], 'year': [k], 'ucode': [l], 'fcasum': [fcasum], 'fcave': [fcave], 'fcstd': [fcstd]}
                fluxdf = pd.DataFrame(data=d)
                # Write to file
                fluxdf.to_csv(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '_fcincorrnoifl.csv')
                # Delete shapefile
                driver = ogr.GetDriverByName("ESRI Shapefile")
                driver.DeleteDataSource(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '.shp')
            else:
                # Delete shapefile
                driver = ogr.GetDriverByName("ESRI Shapefile")
                driver.DeleteDataSource(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '.shp')