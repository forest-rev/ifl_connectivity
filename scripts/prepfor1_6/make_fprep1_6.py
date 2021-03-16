# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 13:30:39 2019
Make scripts for preparing gfc1.6 forest cover for corridor processing.
Run using python command line. E.g.
python ./code/scripts/prepfor1_6/make_fprep1_6.py
@author: pj276
"""

#%% Modules
import fileinput, os, commands
from shutil import copyfile

#%%
# 1.	hf2_reclass_lossyear_submit.sh
# 2.	hf_prepare_forest_submit.sh
# 4.	hf_mosaic_vrt_submit.sh
# gdalbuildvrt for each of the forest layers?
# 8.	hf_water_opening_submit.sh
# 9.	hf2_make_exp_cs_submit.sh
# 10.	hf2_aggregate_rasters_submit.sh

#%%
# Monsoon git dir
gdir = '/home/pj276/projects/ifl_corridors'

# Source script directory
ssdir = gdir + '/code/scripts'

# Basename of parameter file
pbase = 'p_gfc16_111719'

#%%
# Make GFC directories if they don't exist
# GFC directory
gfcdir = '/scratch/pj276/ifl_corridors/gfc/gfc1_6'
# Directories to be created
dlist = ['/gfc1_6_af_ss', '/gfc1_6_af_ss_emt', '/gfc1_6_as_ss', '/gfc1_6_as_ss_emt', '/gfc1_6_sa_ss', '/gfc1_6_sa_ss_emt']
# Loop through and create
for i in dlist:
    if not os.path.exists(gfcdir + i):
        os.makedirs(gfcdir + i)
#%%
# Main directory holding param files and scripts
mdir = gdir + '/code/scripts/prepfor1_6'
# Make dir if it doesn't exist already
if not os.path.exists(mdir):
    os.mkdir(mdir)

#%%
#########################################################################
#---- Section 1 parameter files ----
# Template for parameter files (one for each continent)
ptwrite = """# -*- coding: utf-8 -*- 
\"\"\"Created on Sun Nov 17 13:48:00 2019
This is a parameter file to make it easier to re-run corridor operations.
@author: pj276
\"\"\"

# Libraries
import os, gdal

# Continent code
continent = '{continent}'

# Base directory
tld = '/scratch/pj276/ifl_corridors'
# Tree cover directory
fcd = '/gfc/gfc1_6/gfc1_6_' + continent
# Projected tree data directory
pfdd = '/gfc/gfc1_6/gfc1_6_' + continent +'_ss'
# Edge matched tile directory
emtd = '/gfc/gfc1_6/gfc1_6_' + continent + '_ss_emt'

###############################################################################
# hf2_reclass_lossyear.py
###############################################################################
# Year through which to calculate forest loss or after which to calculate
# forest loss
uyear = '2018'
# Loss year value to use in masking lossyear raster
lyv = 18
# Reclassify only values greater than threshold or less than threshold
# Less than includes upper threshold. Lower threshold is hard coded to zero.
gtorlt = 'lt'
# Reclass loss year name
rclyname = 'loss20002018'
###############################################################################

###############################################################################
#hf_prepare_forest.py
###############################################################################
# Set up arguments
# List of search strings
foresttiles = ['treecover2018']
#foresttiles = ['lybytc2000_' + str(i) + '_' for i in range(0,19,1)]
#['loss20132018', 'loss20002018']
#['datamask', 'gain', 'lossyear', 'treecover2000', 'treecover2012_', 'treecover2012gn']

# Output directory for projected forest data
odfor = tld + pfdd

# File for reference projection
if continent == 'as':
    #as
    prjrast = '/scratch/pj276/ifl_corridors/prjfiles' + '/Hansen_GFC-2016-v1.4_datamask_20N_070E_ss_as.tif'
if continent == 'af':
    #af
    prjrast = '/scratch/pj276/ifl_corridors/prjfiles' + '/Hansen_GFC-2016-v1.4_datamask_10N_020W_ss_af.tif'
if continent == 'sa':
    #sa
    prjrast = '/scratch/pj276/ifl_corridors/prjfiles' + '/Hansen_GFC-2016-v1.4_datamask_30N_100W_ss_sa.tif'

# Projection suffix (e.g. sinusoidal = ss, lambert azimuthal equal area = lazea) and file extension
projsuffix = "_ss.tif"

# Set snap coordinates
xint = 30 # cell size
yint = 30 # cell size

#ss
if continent == 'af':
    # af
    snap_x = -5009377.085697310976684 + (10000*xint)
    snap_y = 1105854.833234372083098 - (10000*yint)
if continent == 'as':
    #as
    snap_x = -4933771.383066884241998 + (10000*xint)
    snap_y = 2212366.254171633161604 - (10000*yint)
if continent == 'sa':
    #sa
    snap_x = -3662648.020845168270171 + (10000*xint)
    snap_y = 3320113.397940382361412 - (10000*yint)

###############################################################################
#hf2_lossyear_by_tc2000.py
###############################################################################
# Individual loss year by tree cover 2000 name
lybytc2000name = 'lybytc2000'
###############################################################################

###############################################################################
#hf_update_forest_cover.py
###############################################################################
# Important arguments for hf_update_forest_cover.py
# Year 2000 tree cover file list
fc2000list = tld + fcd + '/' + 'findextreecov_' + continent + '.txt'
# Gain list
gainlist = tld + fcd + '/' + 'findexgain_' + continent + '.txt'
# Year through which to calculate forest loss
uyear = '2018'
# Loss year value to use in masking lossyear raster
lyv = 12
# If True, use year 2012 gain to update forest cover, if False, just use loss
# through specified year
usegain = False
###############################################################################

###############################################################################
#hf_make_fnet.py
###############################################################################
# Fishnet file name for organinzing files
fnet = tld + '/fishnets' + '/fnet_' + continent + '_ss_v2.shp' # v2 based on gfc15
# Input directory
infnetdir = tld + pfdd
# Output directory (edge matched tiles)
outfnetdir = tld + emtd
# Reference projection file
fnetref = tld + '/prjfiles' + '/ss_' + continent + '.prj'
###############################################################################

###############################################################################
# hf2_mosaic_vrt.py
###############################################################################
# Gdal pixel type. Use gdal.GDT_Int16 for lybtc2000, use gdal.GDT_Byte for tree cover and loss year
gdpt = gdal.GDT_Byte
# Numpy pixel type. Use 'int16' for lybtc2000, use uint8 for tree cover and loss year
nppt = 'uint8'
###############################################################################

###############################################################################
# hf2_aggregate_rasters.py
###############################################################################
# Input dir
aggdir = '/scratch/pj276/ifl_corridors/gfc/gfc1_6/gfc1_6_' + continent + '_ss_emt'
# Output dir
aggodir = '/scratch/pj276/ifl_corridors/gfc/gfc1_6/gfc1_6_' + continent + '_ss_emt' + '_90m' 
# Amount by which to expand bounding box dimensions
exvec = [-90, -90, 90, 90]
# New cell dimension
aggdim = "90"
# No data value
ndval = "255"
# Source file name pattern
srcaggname = ['treecover2018']
#[continent + '_rtile_lybytc2000_' + str(i) + '_' for i in range(0,19,1)]
#[continent + '_rtile_lybytc2000_1_']
# Destination file name pattern
dstaggname = ['treecover2018']
#dstaggname = [continent + '_rtile_lybytc2000_' + str(i) + '_90m_' for i in range(0,19,1)]
#[continent + '_rtile_lybytc2000_1_90m_']

"""

for i in ['af','as','sa']:
    context = {"continent":i}
    with open(mdir + '/' + pbase + '_' + i + '.py','w') as myfile:
        myfile.write(ptwrite.format(**context))

###############################################################################
#---- Section 2 ----
# Scripts
#%%
###############################################################################
#---- hf2_reclass_lossyear_submit_xx.sh ----
# Template for writing submit script
ttwrite = """#!/bin/bash
#SBATCH --job-name=rcly
#SBATCH --chdir=/scratch/pj276/ifl_corridors/output
#SBATCH --output=/scratch/pj276/ifl_corridors/logs/rcly_%A_%a.log
#SBATCH --time=00:05:00
#SBATCH --partition=all
#SBATCH --cpus-per-task=2
#SBATCH --mem 32000
#SBATCH --array={njsfupdate}-{njefupdate}

# The number of jobs is I think equal to the number of lines in the
# file index file

module load anaconda/latest

echo "File ID is:"

echo ${SLURM_ARRAY_TASK_ID}

srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepfor1_6/hf2_reclass_lossyear_{continent}.py ${SLURM_ARRAY_TASK_ID}
"""

# Get number of files in each continent folder to make submit scripts
# and write into the template and then write to file
for k in ['af','as','sa']:
    d = gfcdir + '/' + os.path.basename(gfcdir) + '_' + k
    # List year 2000 tree cover files
    fl = os.listdir(d)
    fl = [i for i in fl if 'lossyear' in i and i.endswith('.tif')]
    # Get number of jobs for reclass loss year script using length of file list
    njsfupdate = 1
    njefupdate = len(fl)
    # Set task ID variable (this is a nuisance, can't figure out how to ignore it in the template)
    SLURM_ARRAY_TASK_ID = 'SLURM_ARRAY_TASK_ID'
    
    context = {
     "njsfupdate": njsfupdate, 
     "njefupdate": njefupdate,
     "SLURM_ARRAY_TASK_ID": SLURM_ARRAY_TASK_ID,
     "continent": k
    }
    with open(mdir + '/hf2_reclass_lossyear_submit_' + k + '.sh','w') as myfile:
        myfile.write(ttwrite.format(**context))

#%%
#---- Make hf2_reclass_lossyear_xx.py scripts ----
for k in ['af','as','sa']:
    copyfile(ssdir + '/hf2_reclass_lossyear.py', mdir + '/hf2_reclass_lossyear_' + k + '.py')
    f = fileinput.FileInput(mdir + '/hf2_reclass_lossyear_' + k + '.py', inplace=True)
    for line in f:
        line = line.replace("import param_file as p1", "import " + pbase + "_" + k + " as p1")
        line = line.replace("sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')",
                            "sys.path.append('/home/pj276/projects/ifl_corridors/code/scripts/prepfor1_6/')")
        print line,
    f.close()

#%%
##---- SUBMIT JOBS ----
#for k in ['af', 'as', 'sa']:
##for k in ['af']:    
#    cmd = "sbatch " +   mdir + '/hf2_reclass_lossyear_submit_' + k + '.sh'
#    print "Submitting Job " + k + " with command: %s" % cmd
#    status, jobnum = commands.getstatusoutput(cmd)
#    if (status == 0):
#        print "Job1 is %s" % jobnum
#    else:
#        print "Error submitting Job " + k

#%%
#########################################################################
#---- hf_update_forest_cover_submit_xx.sh ----
# Template for writing submit script
ttwrite = """#!/bin/bash
#SBATCH --job-name=updfc
#SBATCH --chdir=/scratch/pj276/ifl_corridors/output
#SBATCH --output=/scratch/pj276/ifl_corridors/logs/job_%A_%a.log
#SBATCH --time=00:05:00
#SBATCH --partition=all
#SBATCH --cpus-per-task=2
#SBATCH --mem 32000
#SBATCH --array={njsfupdate}-{njefupdate}

# The number of jobs is I think equal to the number of lines in the
# file index file

module load anaconda/latest

echo "File ID is:"

echo ${SLURM_ARRAY_TASK_ID}

srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepfor1_6/hf_update_forest_cover_{continent}.py ${SLURM_ARRAY_TASK_ID}
"""

# Get number of files in each continent folder to make submit scripts
# and write into the template and then write to file
for k in ['af','as','sa']:
    d = gfcdir + '/' + os.path.basename(gfcdir) + '_' + k
    # List year 2000 tree cover files
    fl = os.listdir(d)
    fl = [i for i in fl if 'treecover2000' in i and i.endswith('.tif')]
    # Get number of jobs for update forest cover script using length of file list
    njsfupdate = 1
    njefupdate = len(fl)
    # Set task ID variable (this is a nuisance, can't figure out how to ignore it in the template)
    SLURM_ARRAY_TASK_ID = 'SLURM_ARRAY_TASK_ID'
    
    context = {
     "njsfupdate": njsfupdate, 
     "njefupdate": njefupdate,
     "SLURM_ARRAY_TASK_ID": SLURM_ARRAY_TASK_ID,
     "continent": k
    }
    with open(mdir + '/hf_update_forest_cover_submit_' + k + '.sh','w') as myfile:
        myfile.write(ttwrite.format(**context))

#%%
##
#---- Make hf_update_forest_cover_xx.py scripts ----
for k in ['af','as','sa']:
    copyfile(ssdir + '/hf_update_forest_cover.py', mdir + '/hf_update_forest_cover_' + k + '.py')
    f = fileinput.FileInput(mdir + '/hf_update_forest_cover_' + k + '.py', inplace=True)
    for line in f:
        line = line.replace("import param_file as p1", "import " + pbase + "_" + k + " as p1")
        line = line.replace("sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')",
                            "sys.path.append('/home/pj276/projects/ifl_corridors/code/scripts/prepfor1_6/')")
        print line,
    f.close()

#%%
#---- SUBMIT JOBS ----
#for k in ['af', 'as', 'sa']:
##for k in ['af']:  
#    cmd = "sbatch " +   mdir + '/hf_update_forest_cover_submit_' + k + '.sh'
#    print "Submitting Job " + k + " with command: %s" % cmd
#    status, jobnum = commands.getstatusoutput(cmd)
#    if (status == 0):
#        print "Job1 is %s" % jobnum
#    else:
#        print "Error submitting Job " + k

#%%
#########################################################################
#---- hf_prepare_forest_submit_xx.sh ----
ttwrite = """#!/bin/bash
#SBATCH --job-name=prep_forest
#SBATCH --chdir=/scratch/pj276/ifl_corridors/output
#SBATCH --output=/scratch/pj276/ifl_corridors/logs/job_%A_%a.log
#SBATCH --time=01:30:00
#SBATCH --partition=all
#SBATCH --cpus-per-task=24
#SBATCH --mem 16000

# should need 2 hours with 12 cpus and mem at 16G
# so maybe 1 hour with 24 cpus
module load anaconda/latest
#module load gdal -this seems to lead to a broken proj4 library

srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepfor1_6/hf_prepare_forest_{continent}.py
"""
# Get number of files in each continent folder to make submit scripts
# and write into the template and then write to file
for k in ['af','as','sa']:
    context = {"continent": k}
    with open(mdir + '/hf_prepare_forest_submit_' + k + '.sh','w') as myfile:
        myfile.write(ttwrite.format(**context))

#%%
#---- Make hf_prepare_forest_xx.py scripts ----
for k in ['af','as','sa']:
    copyfile(ssdir + '/hf_prepare_forest.py', mdir + '/hf_prepare_forest_' + k + '.py')
    f = fileinput.FileInput(mdir + '/hf_prepare_forest_' + k + '.py', inplace=True)
    for line in f:
        line = line.replace("import param_file as p1", "import " + pbase + "_" + k + " as p1")
        line = line.replace("sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')",
                            "sys.path.append('/home/pj276/projects/ifl_corridors/code/scripts/prepfor1_6/')")
        print line,
    f.close()

###---- SUBMIT JOBS ----
#for k in ['af', 'as', 'sa']:
#    cmd = "sbatch " +   mdir + '/hf_prepare_forest_submit_' + k + '.sh'
#    print "Submitting Job " + k + " with command: %s" % cmd
#    status, jobnum = commands.getstatusoutput(cmd)
#    if (status == 0):
#        print "Job1 is %s" % jobnum
#    else:
#        print "Error submitting Job " + k

#%%
#########################################################################
#---- hf2_lossyear_by_tc2000_submit_xx.sh ----
# Template for writing submit script
ttwrite = """#!/bin/bash
#SBATCH --job-name=ly00
#SBATCH --chdir=/scratch/pj276/ifl_corridors/output
#SBATCH --output=/scratch/pj276/ifl_corridors/logs/ly00_%A_%a.log
#SBATCH --time=00:30:00
#SBATCH --partition=all
#SBATCH --cpus-per-task=2
#SBATCH --mem 64000
#SBATCH --array={njsfupdate}-{njefupdate}

module load anaconda/latest

echo "File ID is:"

echo ${SLURM_ARRAY_TASK_ID}

srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepfor1_6/hf2_lossyear_by_tc2000_{continent}.py ${SLURM_ARRAY_TASK_ID}
"""

# Get number of files in each continent folder to make submit scripts
# and write into the template and then write to file
for k in ['af','as','sa']:
    d = gfcdir + '/' + os.path.basename(gfcdir) + '_' + k + '_ss'
    # List year 2000 tree cover files
    fl = os.listdir(d)
    fl = [i for i in fl if 'lossyear' in i and i.endswith('.tif')]
    # Get number of jobs for reclass loss year script using length of file list
    njsfupdate = 1
    njefupdate = len(fl)
    # Set task ID variable (this is a nuisance, can't figure out how to ignore it in the template)
    SLURM_ARRAY_TASK_ID = 'SLURM_ARRAY_TASK_ID'
    
    context = {
     "njsfupdate": njsfupdate, 
     "njefupdate": njefupdate,
     "SLURM_ARRAY_TASK_ID": SLURM_ARRAY_TASK_ID,
     "continent": k
    }
    with open(mdir + '/hf2_lossyear_by_tc2000_submit_' + k + '.sh','w') as myfile:
        myfile.write(ttwrite.format(**context))

#%%
#---- Make hf2_lossyear_by_tc2000_xx.py scripts ----
for k in ['af','as','sa']:
    copyfile(ssdir + '/hf2_lossyear_by_tc2000.py', mdir + '/hf2_lossyear_by_tc2000_' + k + '.py')
    f = fileinput.FileInput(mdir + '/hf2_lossyear_by_tc2000_' + k + '.py', inplace=True)
    for line in f:
        line = line.replace("import param_file as p1", "import " + pbase + "_" + k + " as p1")
        line = line.replace("sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')",
                            "sys.path.append('/home/pj276/projects/ifl_corridors/code/scripts/prepfor1_6/')")
        print line,
    f.close()

#%%
##---- SUBMIT JOBS ----
#for k in ['af', 'as', 'sa']:
##for k in ['af']:    
#    cmd = "sbatch " +   mdir + '/hf2_lossyear_by_tc2000_submit_' + k + '.sh'
#    print "Submitting Job " + k + " with command: %s" % cmd
#    status, jobnum = commands.getstatusoutput(cmd)
#    if (status == 0):
#        print "Job1 is %s" % jobnum
#    else:
#        print "Error submitting Job " + k

#%%
#################################################################################################
#---- hf_mosaic_vrt_submit_xx.sh ----
ttwrite = """#!/bin/bash
#SBATCH --job-name=tile_mos
#SBATCH --chdir=/scratch/pj276/ifl_corridors/output
#SBATCH --output=/scratch/pj276/ifl_corridors/logs/mossa_%A_%a.log
#SBATCH --time=03:30:00
#SBATCH --partition=all
#SBATCH --cpus-per-task=18
#SBATCH --mem 64000

# originally asked for 30 minutes
module load anaconda/latest
module load gdal
srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepfor1_6/hf_mosaic_vrt_{continent}.py
"""
# Write template to file
for k in ['af','as','sa']:
    context = {"continent": k}
    with open(mdir + '/hf_mosaic_vrt_submit_' + k + '.sh', 'w') as myfile:
        myfile.write(ttwrite.format(**context))

#%%
#---- Make hf_mosaic_vrt_xx.py scripts ----
for k in ['af','as','sa']:
    copyfile(ssdir + '/hf_mosaic_vrt.py', mdir + '/hf_mosaic_vrt_' + k + '.py')
    f = fileinput.FileInput(mdir + '/hf_mosaic_vrt_' + k + '.py', inplace=True)
    for line in f:
        line = line.replace("import param_file as p1", "import " + pbase + "_" + k + " as p1")
        line = line.replace("sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')",
                            "sys.path.append('/home/pj276/projects/ifl_corridors/code/scripts/prepfor1_6/')")
        print line,
    f.close()     

##---- SUBMIT JOBS ----
#for k in ['af', 'as', 'sa']:
#    cmd = "sbatch " +   mdir + '/hf_mosaic_vrt_submit_' + k + '.sh'
#    print "Submitting Job " + k + " with command: %s" % cmd
#    status, jobnum = commands.getstatusoutput(cmd)
#    if (status == 0):
#        print "Job1 is %s" % jobnum
#    else:
#        print "Error submitting Job " + k

#%%
#########################################################################
#---- hf2_aggregate_rasters_xx.sh ----
# Template for writing submit script
ttwrite = """#!/bin/bash
#SBATCH --job-name=aggrast
#SBATCH --chdir=/scratch/pj276/ifl_corridors/output
#SBATCH --output=/scratch/pj276/ifl_corridors/logs/aggrast_%A_%a.log
#SBATCH --time=05:40:00
#SBATCH --partition=all
#SBATCH --cpus-per-task=12
#SBATCH --mem 24000
module load anaconda/latest

srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepfor1_6/hf2_aggregate_rasters_{continent}.py
"""
# Write template to file
for i in ['af','as','sa']:
    # Dic that will be expanded in the format command
    context = {"continent": i}
    with open(mdir + '/hf2_aggregate_rasters_submit_' + i + '.sh','w') as myfile:
        myfile.write(ttwrite.format(**context))

#%%
#---- Make hf2_aggregate_rasters_xx.py scripts ----
for k in ['af','as','sa']:
    copyfile(ssdir + '/hf2_aggregate_rasters.py', mdir + '/hf2_aggregate_rasters_' + k + '.py')
    f = fileinput.FileInput(mdir + '/hf2_aggregate_rasters_'+ k + '.py', inplace=True)
    for line in f:
        line = line.replace("import param_file as p1", "import " + pbase + '_' + k + " as p1")
        line = line.replace("sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')",
                            "sys.path.append('/home/pj276/projects/ifl_corridors/code/scripts/prepfor1_6/')")
        print line,
    f.close()

##---- SUBMIT JOBS ----
#for k in ['af', 'as', 'sa']:
#    cmd = "sbatch " +   mdir + '/hf2_aggregate_rasters_submit_' + k + '.sh'
#    print "Submitting Job " + k + " with command: %s" % cmd
#    status, jobnum = commands.getstatusoutput(cmd)
#    if (status == 0):
#        print "Job1 is %s" % jobnum
#    else:
#        print "Error submitting Job " + k
##########################################################################