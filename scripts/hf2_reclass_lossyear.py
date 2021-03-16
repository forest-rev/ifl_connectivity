# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 13:48:48 2019
This script reclassifies lossyear > 12 to 1, 0 otherwise.
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

# Job array index
j = int(sys.argv[1]) - 1 # subtract for python zero indexing

# Loss year value
lyv = p1.lyv #12

# Threshold type
gtorlt = p1.gtorlt

# Reclass loss year name
rclyname = p1.rclyname

# List loss year files
fl = os.listdir(p1.tld + p1.fcd)
fl = [i for i in fl if 'lossyear' in i and i.endswith('.tif')]

# Get file name to work on
fnamely = p1.tld + p1.fcd + '/' + fl[j]

# Loss year
lyear = hfu.raster2array(fnamely)
lya = lyear[0]
del lyear

if gtorlt == 'gt':
    # Reclass loss year > 12 to 1, else 0
    lya = np.where(lya > lyv, 1, 0)
    lya = lya.astype('int8')
else:
    # Reclass loss year <= 12 to 1, else 0
    lya = np.where((lya > 0) & (lya <= lyv), 1, 0)
    lya = lya.astype('int8')

# Output reclassed loss file name
ofn = os.path.dirname(fnamely) + os.sep + os.path.basename(fnamely).replace('lossyear', rclyname)
# Write to array
hfu.array2raster(ofn, fnamely, gdal.GDT_Byte, 'GTiff', lya)
