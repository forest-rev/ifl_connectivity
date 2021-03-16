# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 21:34:56 2019
This script reads in a csv of patch pairs and saves a convex hull shapefile.
@author: pj276
"""

#%%
import pandas as pd
import sys, os, osr, ogr, fiona, subprocess, glob, shutil
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.tools.mask import mask
sys.path.append('/home/pj276/projects/ifl_corridors/code/functions/')
import hf_utilities as hfu
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
#import p1_90m_110718_as as p1
import param_file as p1

#%%
#'/projects/above_gedi/pjantz/ifl_corridors_data/outputs/regular_points/rp10ks1/ptnbs_10k.csv'
tdist = p1.tdist
# Convex hull directory
chulldir = p1.chulldir
#300000 # maximum point distance to consider, in meters
nbcsvall = p1.nbcsvall
# Hinterland forest patches
patches = p1.inpoly
# Directory
odir = p1.odir
# Log file directory
ppc4redodir = p1.ppc4redodir
# Tree cover vrt
tcvrt = p1.tcvrt
# Areas data frame name prefix
adfname_temp = p1.adfname_temp
# Get x and y cell sizes
xint = p1.xint
yint = p1.yint
# Get year
ayear = p1.ayear
# Job array index
j = int(sys.argv[1]) - 1 # subtract for python zero indexing
adfname_temp = adfname_temp + str(j) + '.csv'

#%%
# Skip if output exists
if not os.path.exists(adfname_temp):
    #%%
    # Build a vrt if it doesn't exist
    if not os.path.exists(tcvrt):
        alist = glob.glob(os.path.dirname(tcvrt) + '/' + '*treecover' + str(ayear) + '*')
        with open(os.path.dirname(tcvrt) + '/' + 'tempfilelist_' + str(j) + '.txt', 'a') as tfile:
            for k in alist:
                tfile.write(k + '\n')
        subprocess.call(["gdalbuildvrt", "-input_file_list", os.path.dirname(tcvrt) + '/' + 'tempfilelist_' + str(j) + '.txt', tcvrt])
        os.remove(os.path.dirname(tcvrt) + '/' + 'tempfilelist_' + str(j) + '.txt')
    
    #%%
    # Draw a convex hull around each patch's set of points
    # and calculate subpixel canopy area for each patch by intersecting
    # patches and convex hulls, multiplying canopy cover percent by 
    # 30x30m (or 90x90m) cell area and summing.
    
    # Get projection info from input shapefile
    f = open(os.path.dirname(patches) + os.sep + os.path.basename(patches).replace('shp','prj'))
    fr = f.read() # Read wkt to variable
    f = None # Close file
    srs=osr.SpatialReference(wkt=fr)
    
    # Read in csv
    rdf = pd.read_csv(nbcsvall)
    gppd = rdf.groupby(['gc1','gc2'])
    rdf = None
    gpk = gppd.groups.keys()
    # Subset data frame
    aname = gpk[j]
    group = gppd.get_group(aname)
    
    # Empty lists to hold patch area/quality estimates
    pid1 = [] # patch 1 id
    pid2 = [] # patch 2 id
    
    aname = [str(i) for i in aname]
    # Collapse name to string
    aname = '_'.join(aname)
    # List coordinates as x,y tuples
    tuplist = zip(group.cx1.tolist(), group.cy1.tolist())
    tuplist.extend(zip(group.cx2.tolist(), group.cy2.tolist()))
    # Make dir to hold temporary points, deleting dir first if it exists
    if os.path.exists(odir + chulldir + os.sep + aname):
        shutil.rmtree(odir + chulldir + os.sep + aname)
    os.mkdir(odir + chulldir + os.sep + aname)
    # Convert coords to point shapefile
    # Temporary point shapefile prefix
    tpoints = odir + chulldir + os.sep + aname + '/temp_point_'
    ptsname = tpoints + aname + '.shp'
    hfu.coordinates2point(tuplist,ptsname,srs)
    # Make buffered convex hull for point set
    och1 = odir + chulldir + os.sep + aname + '/' + 'chull_' + aname + '.shp'
    hfu.convexhull(ptsname, och1, buff=1000)
    # Get patch 1 -----------
    g1 = int(aname.split('_')[0])
    pid1.append(g1) # append to list
    # Use fiona and geopandas to dissolve features if there are multiple ones with the same gridcode
    t1 = gpd.read_file(patches) # get projection info
    myproj = t1.crs
    t1 = None
    reader = fiona.open(patches)
    xx = gpd.GeoDataFrame.from_features((x for x in reader if x['properties']['GRIDCODE']==g1))
    if xx.shape[0] > 1:
        xx = xx.dissolve(by='GRIDCODE', aggfunc='sum')
    xx.crs = myproj
    xx['GRIDCODE'] = g1
    # Fix simple self intersections if necessary
    if xx.is_valid.bool() == False:
        xx = xx.buffer(0)
    # Write polygon to shapefile
    dissoshape1 = odir + chulldir + os.sep + aname + '/patch1_temp_' + str(g1) + '_' + str(j) + '.shp'
    xx.to_file(dissoshape1)
    reader.close()
    
    # Get patch 2 -----------
    g2 = int(aname.split('_')[1])
    pid2.append(g2) # append to list
    # Use fiona and geopandas to dissolve features if there are multiple ones with the same gridcode
    t2 = gpd.read_file(patches) # get projection info
    myproj = t2.crs
    t2 = None
    reader = fiona.open(patches)
    xx = gpd.GeoDataFrame.from_features((x for x in reader if x['properties']['GRIDCODE']==g2))
    if xx.shape[0] > 1:
        xx = xx.dissolve(by='GRIDCODE', aggfunc='sum')
    xx.crs = myproj
    xx['GRIDCODE'] = g2
    # Fix simple self intersections if necessary
    if xx.is_valid.bool() == False:
        xx = xx.buffer(0)
    # Write polygon to shapefile
    dissoshape2 = odir + chulldir + os.sep + aname + '/patch2_temp_' + str(g2) + '_' + str(j) + '.shp'
    xx.to_file(dissoshape2)
    reader.close()
    
    # Get geometry of patch 1 intersection and write to file
    oint1 = odir + chulldir + os.sep + aname + '/patch1_chull_temp_' + str(g1) + '_' + str(j) + '.shp'
    intgeo1 = hfu.intersectionGeom(dissoshape1, och1, field1='uid', field2='uid', write='yes', sfname=oint1)
    
    # Get geometry of patch 2 intersection and write to file
    oint2 = odir + chulldir + os.sep + aname + '/patch2_chull_temp_' + str(g2) + '_' + str(j) + '.shp'
    intgeo2 = hfu.intersectionGeom(dissoshape2, och1, field1='uid', field2='uid', write='yes', sfname=oint2)
   
    # Delete shapefiles
    if os.path.exists(ptsname):
        shpDriver = ogr.GetDriverByName('ESRI Shapefile')
        shpDriver.DeleteDataSource(ptsname)
        shpDriver = None
#    if os.path.exists(oint1):
#        shpDriver = ogr.GetDriverByName('ESRI Shapefile')
#        shpDriver.DeleteDataSource(oint1)
#        shpDriver = None
#    if os.path.exists(oint2):
#        shpDriver = ogr.GetDriverByName('ESRI Shapefile')
#        shpDriver.DeleteDataSource(oint2)
#        shpDriver = None
    if os.path.exists(dissoshape1):
        shpDriver = ogr.GetDriverByName('ESRI Shapefile')
        shpDriver.DeleteDataSource(dissoshape1)
        shpDriver = None
    if os.path.exists(dissoshape2):
        shpDriver = ogr.GetDriverByName('ESRI Shapefile')
        shpDriver.DeleteDataSource(dissoshape2)
        shpDriver = None
    if os.path.exists(och1):
        shpDriver = ogr.GetDriverByName('ESRI Shapefile')
        shpDriver.DeleteDataSource(och1)
        shpDriver = None
    
