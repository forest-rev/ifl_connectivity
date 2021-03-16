# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 14:06:54 2019
This script aggregates csvs of distance estimates between points,
labels special cases, and sums overlapping corridors for each patch pair.
It's a bit confusing but this script sums level 1 corridors. Level 0 corridors
were summed using another script. However, the script does  use level 0 corridor
lists to summarize the appropriate cost distance files.
@author: pj276
"""

#%%
import sys, os, osr, gdal
import pandas as pd
#import geopandas as gpd
import numpy as np
sys.path.append('/home/pj276/projects/ifl_corridors/code/functions/')
import hf_utilities as hfu
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
#import p1_20k_100k_052518 as p1
#import p1_10k_100k_060418 as p1
#import p1_rp_300k_060518 as p1
#import p1_10k_300k_060518 as p1
#import p1_90m_110718_as as p1
import param_file as p1

#%%
# Job array index
j = int(sys.argv[1]) - 1 # subtract for python zero indexing
# Main directory
odir = p1.odir
# Corridor count directory
ccdirp = p1.ccdirp
ccdir = odir + ccdirp
# Corridor unique id file
ucpfile = p1.ucpfile

# Csv file prefix
optilep = p1.optilep
#'lcc_ptile_dists_2000_neg1_'
# Year of corridor
ayear = p1.ayear

# Patch dataset
patches = p1.inpoly
#'/projects/above_gedi/pjantz/ifl_corridors_data/hifo/hland_tropics_ssu_sa_subset.shp'
# Output stats csv basename for all stat values per patch pair
statsoutl1 = p1.statsoutl1
# Output stats csv basename for single stat value per patch pair
statsoutl2 = p1.statsoutl2
#'/projects/above_gedi/pjantz/ifl_corridors_data/outputs/regular_points/rp10ks1/corridors/dist_area_stats_'
# Get spatial reference
prjfile = p1.prjfile
#'/projects/above_gedi/pjantz/ifl_corridors_data/prjfiles/ss_sa.prj'
# Area data frame
adfname = p1.adfname
# Get interaction distance
ixdist = p1.ixdist

#%% Get spatial reference
f = open(prjfile)
fr = f.read() # Read wkt to variable
f = None # Close file
srs=osr.SpatialReference(wkt=fr)

#%%
# Read in csv of unique ids. Each row is a two tuple.
ucpfile = pd.read_csv(ucpfile)
ucpfile = ucpfile.reset_index()

# Get target row as a two tuple of integers
ucpair = (int(ucpfile.loc[j,'pp1']),int(ucpfile.loc[j,'pp2']))

# Get patch pair string
s = '_'
pps = s.join([str(ucpair[0]),str(ucpair[1])])
print pps

# List to hold distances, areas, and area weighted fluxes
slist = []

# List corridor folders
dlist = os.listdir(odir + optilep)

# Make empty list to hold target csv files
tflist = [k for k in dlist if pps + '_' in k]
# Grab csv files for the target patch pair
tflist = [pd.read_csv(odir + optilep + '/' + i) for i in tflist]
# Concatenate them
fdf = pd.concat(tflist, ignore_index=True)
# Select only rows with non-zero probability of dispersal
fdf = fdf.loc[fdf['apd'] > 0]
# Save to file
statsoutl1 = statsoutl1 + pps + '.csv'
fdf.to_csv(statsoutl1) 

# If the final output stat file exists, delete it
if os.path.exists(statsoutl2 + pps + '.csv'):
    os.remove(statsoutl2 + pps + '.csv')
    
# Convert probability of disperal list to array
darr = np.array(fdf.apd.tolist())

# Read in raster of probability of dispersal sums
corrtifname = ccdir + '/' + pps + '/corr_' + ayear + '_' + pps + '.tif'
csum = hfu.raster2array(corrtifname)[0]

# Calculate root mean square corridor surface roughness which is 
# the standard deviation of the distribution of corridor probabilities. 
cpsd = np.std(csum[csum > 0].astype("float32"))

# Calculate mean absolute deviation from the average corridor probability.
# Similar to cdsd above but less sensitive to high density values.
cpaad = np.mean(np.abs(csum[csum > 0].astype("float32") - np.mean(csum[csum > 0].astype("float32"))))
  
# Get number of pixels with at least one corridor.
ncells = np.sum(csum > 0)

# Normalize probabilities by their sum
wts = csum/np.sum(csum).astype('float32')

# Read in area csv and get areas (either tree canopy area or patch area)
# Areas are in km2
dat = pd.read_csv(adfname)
tca1 = dat[(dat['pid1'] == int(pps.split('_')[0])) & (dat['pid2'] == int(pps.split('_')[1]))].tcarea1.tolist()[0]
tca2 = dat[(dat['pid1'] == int(pps.split('_')[0])) & (dat['pid2'] == int(pps.split('_')[1]))].tcarea2.tolist()[0]

# Set output name for flux surface
otiff = odir + '/cfinal/flux/corr_' + ayear + '_' + pps + '_flux.tif' # output tiff name

# To visualize use either exponential decay probabilities (25km median dispersal distance)
# and no square root weighting of forest area within patches
# or
# inverse distance times mean of forest area within patches
spatialflux = wts*np.mean(darr).astype('float32')*tca1*tca2
#spatialflux = wts*(1/np.mean(draw).astype('float32'))*tca1*tca2
hfu.array2raster(otiff, corrtifname, gdal.GDT_Float32,'GTiff', spatialflux)

# Set output name for normalized probabilities
ncout = odir + '/cfinal/nprob/corr_' + ayear + '_' + pps + '_nprob.tif' # output tiff name
if not os.path.exists(odir + '/cfinal/nprob'):
    os.mkdir(odir + '/cfinal/nprob')
hfu.array2raster(ncout, corrtifname, gdal.GDT_Float32,'GTiff', wts)
 
# Mean and median probabilities of dispersal
aveprob = np.mean(darr).astype('float32')
medprob = np.median(darr).astype('float32')

# Mean and median flux
aveflux = aveprob*tca1*tca2
medflux = medprob*tca1*tca2

# Inverse distance times mean of areas
#flux2 = (1/np.average(draw).astype('float32'))*np.sqrt((pahuArea1/1000000))*np.sqrt((pahuArea2/1000000))
#flux2 = (1/np.average(draw).astype('float32'))*np.sqrt((tca1))*np.sqrt((tca2))
#flux2 = (1/np.average(draw).astype('float32'))*tca1*tca2

# Inverse distance squared times square roots of areas
#flux3 = (1/np.average(draw).astype('float32')**2)*np.sqrt((pahuArea1/1000000))*np.sqrt((pahuArea2/1000000))
#flux3 = (1/np.average(draw).astype('float32')**2)*np.sqrt((tca1))*np.sqrt((tca2))
#flux3 = (1/np.average(draw).astype('float32')**2)*tca1*tca2

slist.append([int(pps.split('_')[0]),int(pps.split('_')[1]),aveprob,medprob,tca1,tca2,aveflux,medflux,cpsd,cpaad,ncells])

# Write stats to dataframe    
# Dump to data frame
nbdf = pd.DataFrame(slist)
# Name columns
nbdf.columns = ['gc1', 'gc2', 'p10probave', 'p10probmed', 'gc1Area', 'gc2Area', 'p10fluxave', 'p10fluxmed', 'cprobVar1', 'cprobVar2', 'ncells']
# Name index
nbdf.index.name = 'rowid'
# Write to file
statsoutl2 = statsoutl2 + pps + '.csv'
nbdf.to_csv(statsoutl2)    


