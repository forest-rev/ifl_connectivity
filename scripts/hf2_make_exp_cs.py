# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 18:24:03 2018
This script converts forest cover to a cost surface using an exponential
transformation from Keeley et al. 2016.
@author: pj276
"""
# Import packages
import numpy as np
import sys, os, gdal
sys.path.append('/home/pj276/projects/ifl_corridors/code/functions/')
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
# For hf_utilities, the modules are available within the functions
# when they're used. Otherwise, one needs to prefix with hfu.
# e.g. hfu.shp.Writer
import hf_utilities as hfu
import param_file as p1

#----
# Job array index.
j = int(sys.argv[1])

#----
# Input directory
idir = p1.expidir #'/scratch/pj276/gfc/gfc1_4_as_ss_emt'
# Output directory
odir = p1.expodir #'/scratch/pj276/gfc/gfc1_4_as_ss_cs'
# Morphology directory
mdir = p1.expmdir #'/scratch/pj276/gfc/gfc1_4_morph'
# Morphology neighborhood (in # of cells)
mnbo = p1.mnbo #'33'
# Region
reg = p1.reg
# Shape parameter for exponential function.
# Use -0.01 to approximate a linear transformation from forest cover to suitability.
# Use -2 to approximate a non-linear transformation where suitability drops off steeply with forest cover decrease
# Use 2 to approximate a non-linear transformation where suitability drops off gradually with forest cover decrease
cshape = p1.cshape #-0.01
# Name corresponding to shape parameter
spname = p1.spname #'neg01' #'neg1'
# Year
y = p1.expy #2012
# Switch for using gain or not
usegain = p1.usegain
# Need to do this to make sure proper suffix is added after year in file name for later year tree cover
if usegain == 'yes' and y == 2012:
    fpat = 'gn_'
else:
    fpat = '_'
# Switch for using morphology or not
usemorph = p1.usemorph #'no'

#----
if usemorph == 'yes':
    # Check for forest cover tile and make sure output doesn't exist
    if (os.path.exists(idir + '/' + reg + '_rtile_treecover' + str(y) + fpat + str(j) + '.tif')) and not os.path.exists(odir + '/' + reg + '_rtile_cs' + str(y) + spname + '_' + str(j) + '.tif'):
        # Get forest cover tile
        fortile = idir + '/' + reg + '_rtile_treecover' + str(y) + fpat + str(j) + '.tif'
        fco = hfu.raster2array(fortile)
        fc = fco[0]
        # Burn in MSPA classes
        mspatile = mdir + '/mspa_' + str(y) + '_' + mnbo + '_' + str(j) + '.tif'
        mspa = hfu.raster2array(mspatile)[0]
        # Codes: core-17; perforation-5,37,69; edge-3,35,67; bridge-33; loop-65; branch-1, patch-9
        minusarray = np.zeros(mspa.shape)
        minusarray[(mspa == 5) | (mspa == 37) | (mspa == 69)] = 5
        minusarray[(mspa == 3) | (mspa == 35) | (mspa == 67)] = 10
        minusarray[(mspa == 33) | (mspa == 65) | (mspa == 1)] = 15
        minusarray[mspa == 9] = 20
        fc = fc - minusarray
        fc[fc < 0] == 0 # Recode any negative values
        mspa = None
        minusarray = None
        # Transform forest cover values (high forest cover values = low cost)
        fc = fc.astype(np.float32)/100.0 # Convert to 0-1 range
        # Use exponential tranformation from Keeley et al. 2016
        # Transformed values range from 1 (low resistance) to 100 (high resistance)
        fc = 100.0-99.0*((1-np.exp(-cshape*fc))/(1-np.exp(-cshape)))
        # Get large water body tile
        rdist = 'r16' # radius of structuring element, in cells, used to create large water body mask
        lwbtile = idir + '/' + reg + '_water_open_' + rdist + '_' + str(j) + '.tif'
        # Read in large water body tile to array
        lwb = hfu.raster2array(lwbtile)[0]
        # Burn in 255 for large water bodies.
        fc[lwb==1] = 255
        lwb = None
        # Burn in no data value as 255 from data mask
        dmasktile = idir + '/' + reg + '_rtile_datamask_' + str(j) + '.tif'
        dmask = hfu.raster2array(dmasktile)[0]
        fc[dmask==0] = 255
        dmask = None
        # Write to file
        csrast = odir + '/' + reg + '_rtile_cs' + str(y) + spname + '_' + str(j) + '.tif'
        # Save as byte
        hfu.array2raster(csrast,fortile,gdal.GDT_Byte,'GTiff',fc)
else:
    # Check for forest cover tile and make sure output doesn't exist
    if (os.path.exists(idir + '/' + reg + '_rtile_treecover' + str(y) + fpat + str(j) + '.tif')) and not os.path.exists(odir + '/' + reg + '_rtile_cs' + str(y) + spname + '_' + str(j) + '.tif'):
        # Get forest cover tile
        fortile = idir + '/' + reg + '_rtile_treecover' + str(y) + fpat + str(j) + '.tif'
        fco = hfu.raster2array(fortile)
        fc = fco[0]
        # Transform forest cover values (high forest cover values = low cost)
        fc = fc.astype(np.float32)/100 # Convert to 0-1 range
        # Use exponential tranformation from Keeley et al. 2016
        # Transformed values range from 1 (low resistance) to 100 (high resistance)
        fc = 100.0-99.0*((1-np.exp(-cshape*fc))/(1-np.exp(-cshape)))
        # Get large water body tile
        rdist = 'r16' # radius of structuring element, in cells, used to create large water body mask
        lwbtile = idir + '/' + reg + '_water_open_' + rdist + '_' + str(j) + '.tif'
        # Read in large water body tile to array
        lwb = hfu.raster2array(lwbtile)[0]
        # Burn in 255 for large water bodies.
        fc[lwb==1] = 255
        lwb = None
        # Burn in no data value as 255 from data mask
        dmasktile = idir + '/' + reg + '_rtile_datamask_' + str(j) + '.tif'
        dmask = hfu.raster2array(dmasktile)[0]
        fc[dmask==0] = 255
        dmask = None
        # Write to file
        csrast = odir + '/' + reg + '_rtile_cs' + str(y) + spname + '_' + str(j) + '.tif'
        # Save as byte
        hfu.array2raster(csrast,fortile,gdal.GDT_Byte,'GTiff',fc)
