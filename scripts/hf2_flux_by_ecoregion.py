# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 23:02:14 2019
Calculate flux by ecoregion
@author: pj276
"""

#%%
# Imports
import subprocess, geopandas as gpd, rasterio, os, ogr, sys
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

# Directory holding original flux files
ofluxdir = p1.ofluxdir

# Directory that will hold output flux files
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

# Specify directory holding zone shapefile
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

# Loop through and calculate flux by zone
for k in year:
    # Get tile index
    atx = idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/flux/' + j + '_' + i + '_' + k + '_cflux.shp'
    df = gpd.read_file(atx)
    # Get tile index bounds
    minx, miny, maxx, maxy = df.geometry.total_bounds
    b = [box(minx, miny, maxx, maxy)]
    d = {'id': [1], 'region': [i]}
    dd = pd.DataFrame(data=d)
    bbox = gpd.GeoDataFrame(dd, geometry=b)
    bbox.crs = df.crs
    bboxgeog = bbox.to_crs(shpfile.crs)            
    # Loop through unique names/codes
    # Test if ucode is string, treat differently if so
    if isinstance(ucode[jai], str):
        l = ucode[jai]
        lorig = l
        if ' ' in l:
            l = l.replace(' ', '_')
    else:
        l = ucode[jai]
        lorig = l
        l = np.str(np.int(l))
    # Output tiff name
    otiff1 = idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '_csum.tif'
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
        # Check if polygon intersects with current tile index
        minx, miny, maxx, maxy = cx.geometry.total_bounds
        b = [box(minx, miny, maxx, maxy)]
        d = {'id': [1], 'ucode': [l]}
        dd = pd.DataFrame(data=d)
        bboxcx = gpd.GeoDataFrame(dd, geometry=b)
        bboxcx.crs = cx.crs
        ix = gpd.overlay(bboxgeog, bboxcx, how='intersection', use_sindex=True)
        if ix.shape[0] >= 1:
            # Save file if needed
            cx = cx.to_crs(df.crs)
            cx.to_file(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '.shp')              
            # List polygon and tile index intersections using ogr
            driver = ogr.GetDriverByName("ESRI Shapefile")
            file1 = driver.Open(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '.shp',0)
            layer1 = file1.GetLayer()
            feat1 = layer1.GetFeature(0)
            geom1 = feat1.GetGeometryRef()
            hitlist = []
            for z in range(0,df.shape[0]):
                file2 = driver.Open(atx,0)
                layer2 = file2.GetLayer()
                layer2.SetAttributeFilter("FID = " + str(z))
                for q in range( 0, layer2.GetFeatureCount() ):
                    feat2 = layer2.GetNextFeature()
                    geom2 = feat2.GetGeometryRef()
                    if geom1.Intersects(geom2) == 1:
                        print 'item ' + str(z) + ' is hit'
                        hitlist.append(feat2.GetFieldAsString(0))
                    else:
                        print 'item ' + str(z) + ' is no hit'
            # Process further if there are intersections
            if len(hitlist) > 0:
                # Save text file of intersecting flux file names            
                iflist = idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '_cflux.txt'
                # Path to cflux
                path = idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + ofluxdir + '/'
                # Concat path and file name
                hitlist = [path + n for n in hitlist]
                # Write to file
                with open(iflist, 'w') as text_file:
                    for item in hitlist:
                        text_file.write('%s\n' % item)
                # Build vrt from text file
                avrt = idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '_cflux.vrt'
                subprocess.call(["gdalbuildvrt", "-vrtnodata", "0", "-separate", "-input_file_list", iflist, avrt])
                # Sum flux files intersecting with polygon boundary
                csum = hfu.rastersum(avrt)
                # Write to file
                hfu.array2raster(otiff1, avrt, gdal.GDT_Float32, 'GTiff', csum)
                # Get geometry from zone shapefile
                with fiona.open(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '.shp', "r") as shapefile:
                    geoms = [feature["geometry"] for feature in shapefile]
                # Use zone boundary to mask flux sum layer
                with rasterio.open(otiff1) as src:
                    out_image, out_transform = mask(src, geoms, crop=True)
                    out_meta = src.meta.copy()
                # Sum flux values to get total for the zone
                tcsum = np.sum(out_image)
                # Delete rectangle cropped raster
                os.remove(otiff1)
                # Save zone cropped raster
                out_meta.update({"driver": "GTiff",
                                 "compress": "lzw",
                                 "dtype": rasterio.float32,
                                 "height": out_image.shape[1],
                                 "width": out_image.shape[2],
                                 "transform": out_transform})
                with rasterio.open(otiff1, "w", **out_meta) as dest:
                    dest.write(out_image.astype(rasterio.float32))
                # One line dataframe to hold flux and other info
                d = {'region': [i], 'scenario': [j], 'year': [k], 'ucode': [l], 'flux': [tcsum]}
                fluxdf = pd.DataFrame(data=d)
                # Write to file
                fluxdf.to_csv(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '_flux.csv')
