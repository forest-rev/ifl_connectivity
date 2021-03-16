# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 12:16:53 2017
This script uses year forest loss information to update year 2000 canopy
cover. This does not deal with forest growth (e.g. increase in canopy
density).
@author: pj276
"""

# Import packages
import sys, os, gdal
import numpy as np
sys.path.append('/home/pj276/projects/ifl_corridors/code/functions/')
# For hf_utilities, the modules are available within the functions
# when they're used. Otherwise, one needs to prefix with hfu.
# e.g. hfu.shp.Writer
import hf_utilities as hfu
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
sys.path.append('/home/pj276/projects/undp_connectivity/code/parameter_files/')
import param_file as p1

# Use year 2012 gain switch
usegain = p1.usegain

# Job array index
j = int(sys.argv[1]) - 1 # subtract for python zero indexing

# Update year
uyear = p1.uyear #'2012'

# Loss year value
lyv = p1.lyv #12

# List year 2000 tree cover files
fl = os.listdir(p1.tld + p1.fcd)
fl = [i for i in fl if 'treecover2000' in i and i.endswith('.tif')]

# Get file name to work on
fname2000 = p1.tld + p1.fcd + '/' + fl[j]

# Use to build file names for forest operations and read to array
# Tree cover in 2000
treecov = hfu.raster2array(fname2000)
tca = treecov[0]

# Build loss year file name from year 2000 forest cover name
fnameloss = os.path.dirname(fname2000) + '/' + os.path.basename(fname2000).replace('treecover2000','lossyear')

# Loss year
lyear = hfu.raster2array(fnameloss)
lya = lyear[0]
del lyear

if usegain == True:
    # Build gain file name from year 2000 forest cover name
    fnamegain = os.path.dirname(fname2000) + '/' + os.path.basename(fname2000).replace('treecover2000','gain')
    
    # Gain
    gain = hfu.raster2array(fnamegain)
    gn = gain[0]
    del gain
    
    # Calculate forest loss through whatever year and update canopy cover.
    treecoverupdate = np.where((lya > 0) & (lya <= lyv), 0, tca)
    # Update with gain using 50% value
    treecoverupdategain = np.where(gn==1, 50, treecoverupdate)
    
    # Make sure int8 data type
    treecoverupdate = treecoverupdate.astype('int8')
    treecoverupdategain = treecoverupdategain.astype('int8')
    
    # Output treecover file name
    ofn = os.path.dirname(fname2000) + os.sep + os.path.basename(fname2000).replace('treecover2000','treecover' + uyear)
    # Write to array
    hfu.array2raster(ofn, fname2000, gdal.GDT_Byte, 'GTiff', treecoverupdate)
    del treecoverupdate
    
    # Output treecover with gain file name
    ofn = os.path.dirname(fname2000) + os.sep + os.path.basename(fname2000).replace('treecover2000','treecover' + uyear + 'gn')
    # Write to array
    hfu.array2raster(ofn, fname2000, gdal.GDT_Byte, 'GTiff', treecoverupdategain)
else:
    
    # Calculate forest loss through whatever year and update canopy cover.
    treecoverupdate = np.where((lya > 0) & (lya <= lyv), 0, tca)
    
    # Make sure int8 data type
    treecoverupdate = treecoverupdate.astype('int8')
    
    # Output treecover file name
    ofn = os.path.dirname(fname2000) + os.sep + os.path.basename(fname2000).replace('treecover2000','treecover' + uyear)
    # Write to array
    hfu.array2raster(ofn, fname2000, gdal.GDT_Byte, 'GTiff', treecoverupdate)
    del treecoverupdate
