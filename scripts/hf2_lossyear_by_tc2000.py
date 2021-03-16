# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 16:23:46 2020
This script steps through the lossyear layer, overlaying each year
with the year 2000 tree cover layer. 
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

# Loss year by tree cover 2000 name
lybytc2000name = p1.lybytc2000name

# List loss year files
fl = os.listdir(p1.tld + p1.pfdd)
fl = [i for i in fl if 'lossyear' in i and i.endswith('.tif')]
fl.sort()
# Get file name to work on
fnamely = p1.tld + p1.pfdd + '/' + fl[j]

# Loss year
lyear = hfu.raster2array(fnamely)
lya = lyear[0]
del lyear

# List tree cover files
y2k = os.listdir(p1.tld + p1.pfdd)
y2k = [i for i in y2k if 'treecover2000' in i and i.endswith('.tif')]
y2k.sort()
# Get file name to work on
y2knamely = p1.tld + p1.pfdd + '/' + y2k[j]
# Year 2000 tree cover
y2k = hfu.raster2array(y2knamely)[0]

# Loop through years
for lyv in range(0, 19, 1):
    # Reclass loss year
    ly = np.where(lya == lyv, 1, 0)
    # Multiply loss year 1/0 raster by percent tree cover by 900m2
    ly = ly*y2k.astype('float32')*900/100.0
    # Output tree cover loss name
    newfname = lybytc2000name + '_' + str(lyv)
    ofn = os.path.dirname(fnamely) + os.sep + os.path.basename(fnamely).replace('lossyear', newfname)
    # Write to array
    hfu.array2raster(ofn, fnamely, gdal.GDT_Int16, 'GTiff', np.rint(ly).astype('int16'))
