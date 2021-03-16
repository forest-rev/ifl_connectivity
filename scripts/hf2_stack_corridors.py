# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 14:15:01 2018
This script stacks corridors.
@author: pj276
"""

# Import packages
import gdal
import subprocess, sys, os, shapefile as shp
import glob
# Append path to functions
sys.path.append('/home/pj276/projects/ifl_corridors/code/functions/')
# For hf_utilities, the modules are available within the functions
# when they're used. Otherwise, one needs to prefix with hfu.
# e.g. hfu.shp.Writer
import hf_utilities as hfu
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
#import p1_10k_300k_060518 as p1
#import p1_90m_110718_as as p1
#import p1_90m_ifl_2013_sa as p1
import param_file as p1

# Working directory
odir = p1.odir
### List of corridor rasters.
##flist = p1.sdir + '/clist.txt'
#flist = '/projects/above_gedi/pjantz/ifl_corridors_data/outputs/regular_points/rp10k300ks1/corridors/clist.txt'
flist = glob.glob(odir + '/cfinal/flux/corr_*')
flist = [x for x in flist if x.endswith('.tif')]

# Job array index. Use this to select a feature from the polygon
j = int(sys.argv[1])-1
# Original index for numbering
k = int(sys.argv[1])

# Get polygon coordinate mins and maxes
fnet = p1.fnet
sf = shp.Reader(fnet)
ashape = sf.shapes()[j]
xpts = [] # list to hold x coords
ypts = [] # list to hold y coords
# Populate lists of x and y coords
for u in ashape.points:
    xpts.append(u[0])
    ypts.append(u[1])

# Get mins and maxes
xmin = str(min(xpts))
xmax = str(max(xpts))
ymin = str(min(ypts))
ymax = str(max(ypts))
       
#------------------------------
# Write file list to text file
with open(odir + '/cfinal/flux/vrtlist_' + str(k) + '.txt', 'w') as text_file:
    for p in flist:
        text_file.write(p + '\n')
# Set output vrt name
vrtcorr = odir + '/cfinal/flux/vrtcorr_' + str(k) + ".vrt"
if os.path.exists(vrtcorr):
    os.remove(vrtcorr)
# Build vrt for tile area
subprocess.call(["gdalbuildvrt", "-separate", "-te", xmin, ymin, xmax, ymax, "-input_file_list", odir + '/cfinal/flux/vrtlist_' + str(k) + '.txt', vrtcorr])
#corrarray = hfu.raster2array3d(vrtcorr)[0]

# Check for no data values
#nantest = np.isnan(corrarray).any()
#corraray = None

# If no no data values, calculate stats
#if nantest == False:
newRasterfn = odir + '/cfinal/flux/sumcorr_' + str(k) + '.tif'
#zz = hfu.tile_mosaic(fnet, vrtcorr, newRasterfn, j, 'sum')
zz = hfu.tile_shingle(vrtcorr, newRasterfn, gdal.GDT_Float32, 'float32', 'sum')

# Clean up
os.remove(vrtcorr)
os.remove(odir + '/cfinal/flux/vrtlist_' + str(k) + '.txt')
