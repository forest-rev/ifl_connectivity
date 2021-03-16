# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 12:52:57 2019
This script loops through convex hull shapefiles representing the area of each patch
involved in connectivity with an adjacent patch. The union of these patches gives
an indication of the total patch area involved in connectivity.
@author: pj276
"""

#%%
# Libraries
import sys, os, glob, ogr, subprocess
import pandas as pd, geopandas as gpd, numpy as np
import fiona, rasterio
from rasterio.tools.mask import mask
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
import param_file as p1
#import p1_90m_ifl_2013_af as p1

#%%
# Directories, files, etc.
odir = p1.odir # main output dir
chulldir = p1.chulldir # convex hull directory
chfname = odir + chulldir + '/convex_hull_area_' # Basename of output for saving convex hull areas

# Tree cover vrt
tcvrt = p1.tcvrt
# Year
ayear = p1.ayear
# Dimensions
xint = p1.xint
yint = p1.yint

# Job array index
j = int(sys.argv[1]) - 1 # subtract for python zero indexing

# Get list of unique patch ids from convex hull shapefiles
fnlist = []
fullnamelist = []
if os.path.exists(odir + chulldir):
    for dir, subdir, files in os.walk(odir + chulldir):
        for fname in files:
            if fname.endswith(".shp") and 'ushape' not in fname:
                fnlist.extend([fname.split('_')[3]])
                fullnamelist.extend([os.path.join(dir + "/", fname)])
# Get unique patch ids and sort
fnlist = list(set(fnlist))
fnlist.sort()

patchid = fnlist[j]
chfname = chfname + patchid + '.csv'

# Get unique full file names
fullnamelist = list(set(fullnamelist))

#%%
# If convex hull and canopy covered areas haven't been calculated at the patch level, run code
if not os.path.exists(chfname):
    # Build a vrt of canopy cover if it doesn't exist
    if not os.path.exists(tcvrt):
        alist = glob.glob(os.path.dirname(tcvrt) + '/' + '*treecover' + str(ayear) + '*')
        with open(os.path.dirname(tcvrt) + '/' + 'tempfilelist.txt', 'a') as tfile:
            for k in alist:
                tfile.write(k + '\n')
        subprocess.call(["gdalbuildvrt", "-input_file_list", os.path.dirname(tcvrt) + '/' + 'tempfilelist.txt', tcvrt])
        os.remove(os.path.dirname(tcvrt) + '/' + 'tempfilelist.txt')
    
# Get shapefiles
shapefiles = [i for i in fullnamelist if '_chull_temp_' + patchid + '_' in i]

# Output lists
pids = list() # empty list to hold patch ids
alist = list() # empty list to hold polygon areas
tcarea = list() # empty list to hold canopy cover areas

#%%
# Get area of unioned shapefiles
# Add shapefiles all to a geodataframe
gdf = pd.concat([gpd.read_file(shp) for shp in shapefiles]).pipe(gpd.GeoDataFrame)
# Union them
su = gdf.unary_union
# Convert geometry to a geodataframe
su = gpd.GeoDataFrame(gpd.GeoSeries(su))
# Set geometry column
su = su.rename(columns={0:'geometry'}).set_geometry('geometry')
# Set projection using first shapfile in the list
su.crs = gpd.read_file(shapefiles[0]).crs
# Get area (convert to km2) from geodataframe                     
a1 = su.area/1000000
a1 = a1.item() # get just the value

# Append id and area to lists
pids.append(patchid)
alist.append(a1)

# Write unioned polygons to file
oufname = odir + chulldir + '/ushape_' + patchid + '.shp'
su.to_file(oufname)

# Use the unioned area to mask the canopy cover layer and extract canopy covered area
# Get polygon masked array for patch 1
with fiona.open(oufname, "r") as shapefile:
    geoms = [feature["geometry"] for feature in shapefile]
with rasterio.open(tcvrt) as src:
    out_image, out_transform = mask(src, geoms, crop=True)
    out_meta = src.meta.copy()
tcsum = np.sum((out_image/100.0)*xint*yint/1000000)
tcarea.append(tcsum)

# Delete unioned shapefile
if os.path.exists(oufname):
    shpDriver = ogr.GetDriverByName('ESRI Shapefile')
    shpDriver.DeleteDataSource(oufname)
    shpDriver = None

#%%
# Save patch ids, unioned areas, and canopy cover area to csv
areasdf = pd.DataFrame({'pid': pids, 'charea': alist, 'tcarea': tcarea})
areasdf = areasdf[['pid','charea','tcarea']]
areasdf.to_csv(chfname,index_label='rownum')
