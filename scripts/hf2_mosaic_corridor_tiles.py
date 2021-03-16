# -*- coding: utf-8 -*-
"""
Created on Mon Jun 04 13:29:10 2018
Mosaic edge matched corridors.
@author: pj276
"""
# Modules
import sys, os, glob, subprocess
sys.path.append('/home/pj276/projects/ifl_corridors/code/functions/')
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
#import p1_10k_300k_060518 as p1
#import p1_90m_110718_as as p1
#import p1_90m_ifl_2000_sa as p1
#import p1_90m_ifl_2013_sa as p1
import param_file as p1

# Get name of working folder
wdir = p1.wdir
# List tiled corridor files and write to file.
txtout = p1.odir + '/ctilelist.txt'
#txtout = '/projects/above_gedi/pjantz/ifl_corridors_data/outputs/regular_points/rp10k300ks1/corridors/ctilelist.txt'
alist = glob.glob(p1.odir + '/cfinal/flux/sumcorr_*')
alist = [x for x in alist if x.endswith('.tif')]
#alist = glob.glob('/projects/above_gedi/pjantz/ifl_corridors_data/outputs/regular_points/rp10k300ks1/corridors/sumcorr_*')
with open(txtout, "w") as text_file:
    for item in alist:
        text_file.write("%s\n" % item)

# Build vrt for tile area
##vrtcorr = p1.idir + '/ctile.vrt'
vrtcorr = p1.odir + '/cfinal/flux/ctile.vrt'
#vrtcorr = '/projects/above_gedi/pjantz/ifl_corridors_data/outputs/regular_points/rp10k300ks1/corridors/ctile.vrt'
subprocess.call(["gdalbuildvrt", "-input_file_list", txtout, vrtcorr])

# Convert vrt to tiff
##otiff = p1.idir + '/csum_' + wdir + '_mos.tif'
otiff = p1.odir + '/cfinal/flux/csum_' + wdir + '_mos.tif'
#otiff = '/projects/above_gedi/pjantz/ifl_corridors_data/outputs/regular_points/rp10k300ks1/corridors/csum_10k300ks1_mos.tif'
subprocess.call(["gdal_translate", "-co", 'COMPRESS=LZW', vrtcorr, otiff])

# Make an integer version
#fras = hfu.raster2array(otiff)[0]
#fras = fras*1000
#fras = fras.astype('int32')
#oint = os.path.dirname(otiff) + '/' + os.path.basename(otiff).split('.')[0] + '_int32.tif'
#hfu.array2raster(oint, otiff, gdal.GDT_Int32, 'GTiff', fras)

# Clean up
#for g in alist:
#    os.remove(g)
os.remove(txtout)
os.remove(vrtcorr)