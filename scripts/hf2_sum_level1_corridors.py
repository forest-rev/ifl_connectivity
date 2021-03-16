# -*- coding: utf-8 -*-
"""
Created on Thu Jul 05 15:00:38 2018
This script sums bundles of rasters up to a first level.
@author: pj276
"""
# Import modules
import sys, os, subprocess, gdal
import pandas as pd

# Append path to functions
sys.path.append('/home/pj276/projects/ifl_corridors/code/functions/')
# For hf_utilities, the modules are available within the functions
# when they're used. Otherwise, one needs to prefix with hfu.
# e.g. hfu.shp.Writer
import hf_utilities as hfu
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
#import p1_90m_110718_as as p1
import param_file as p1

# Set variables
ayear = p1.ayear
odir = p1.odir
sdir = p1.sdir
ucpfile = p1.ucpfile
ccdirp = p1.ccdirp

# Read in csv of 
# Job array index
j = int(sys.argv[1]) - 1 # subtract for python zero indexing

# Read in csv of unique ids. Each row is a two tuple.
ucpfile = pd.read_csv(ucpfile)
ucpfile = ucpfile.reset_index()

# Get target row as a two tuple of integers
ucpair = (int(ucpfile.loc[j,'pp1']),int(ucpfile.loc[j,'pp2']))

# List corridor folders
dlist = os.listdir(sdir)

# Make empty list to hold target files
tflist = []

# Loop through dir list and append tif files if they match the target row ids
for j in dlist:
    if int(j.split('_')[0]) == ucpair[0] and int(j.split('_')[1]) == ucpair[1]:
        myfile = os.listdir(sdir + '/' + j)
        myfile = [i for i in myfile if i.endswith('.tif')]
        if len(myfile) > 0:
            tflist.append(sdir + '/' + j + '/' + myfile[0])

# Build vrt of corridors
corrvrt = odir + '/corr_' + ayear + '_' + str(ucpair[0]) + '_' + str(ucpair[1]) + '.vrt'
# Build a vrt from a file list
acmd = ["gdalbuildvrt", "-vrtnodata", "0", "-separate", corrvrt]
for myfile in tflist:
    acmd.append(myfile)
subprocess.call(acmd)

# Output raster count name. Make dir to hold it.
if not os.path.exists(odir + ccdirp + '/' + str(ucpair[0]) + '_' + str(ucpair[1])):
    os.mkdir(odir + ccdirp + '/' + str(ucpair[0]) + '_' + str(ucpair[1]))
otiff1 = odir + ccdirp + '/' + str(ucpair[0]) + '_' + str(ucpair[1]) + '/corr_' + ayear + '_' + str(ucpair[0]) + '_' + str(ucpair[1]) + '.tif'

# Sum rasters
csum = hfu.rastersum(corrvrt)

# Write to tiff
hfu.array2raster(otiff1, corrvrt, gdal.GDT_Float32, 'GTiff', csum)

# Clean up
os.remove(corrvrt)
