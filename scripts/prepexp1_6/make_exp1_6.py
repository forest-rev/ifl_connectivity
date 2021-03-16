# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 09:55:54 2020
Make scripts for second section (cost surfaces, aggregation to 90m)
@author: pj276
"""

#%% Modules
import commands, fileinput, os
from shutil import copyfile
import numpy as np

#%%
# Monsoon git dir
gdir = '/home/pj276/projects/ifl_corridors'

# Source script directory
ssdir = gdir + '/code/scripts'

# Basename of parameter file
pbase = 'p1_021220_expcs'

#%%
# Make cost surface directories if they don't exist
# GFC directory
gfcdir = '/scratch/pj276/ifl_corridors/gfc/gfc1_6'
# Directories to be created
dlist = ['/gfc1_6_af_ss_cs_neg01', '/gfc1_6_af_ss_cs_neg2', '/gfc1_6_af_ss_cs_pos2',
         '/gfc1_6_as_ss_cs_neg01', '/gfc1_6_as_ss_cs_neg2', '/gfc1_6_as_ss_cs_pos2',
         '/gfc1_6_sa_ss_cs_neg01', '/gfc1_6_sa_ss_cs_neg2', '/gfc1_6_sa_ss_cs_pos2',
         '/gfc1_6_af_ss_cs_neg01_90m', '/gfc1_6_af_ss_cs_neg2_90m', '/gfc1_6_af_ss_cs_pos2_90m',
         '/gfc1_6_as_ss_cs_neg01_90m', '/gfc1_6_as_ss_cs_neg2_90m', '/gfc1_6_as_ss_cs_pos2_90m',
         '/gfc1_6_sa_ss_cs_neg01_90m', '/gfc1_6_sa_ss_cs_neg2_90m', '/gfc1_6_sa_ss_cs_pos2_90m']
# Loop through and create
for i in dlist:
    if not os.path.exists(gfcdir + i):
        os.makedirs(gfcdir + i)
#%%
# Main directory holding param files and scripts
mdir = gdir + '/code/scripts/prepexp1_6'
# Make dir if it doesn't exist already
if not os.path.exists(mdir):
    os.mkdir(mdir)

#%%
#########################################################################
#---- Section 1 parameter file ----
# Template for parameter files (one for each continent and exp scenario and year)
ptwrite = """# -*- coding: utf-8 -*-
\"\"\"Created on Wed Feb 12 09:59:00 2020
This is a parameter file to make it easier to re-run corridor operations.
@author: pj276
\"\"\"

# Libraries
import os

# Continent code
continent = '{continent}'
# Scenario
scen = '{scen}'
# Year
expyear = {expyear}

###############################################################################
# hf2_make_exp_cs.py
###############################################################################
# Input directory
expidir = '/scratch/pj276/ifl_corridors/gfc/gfc1_6/gfc1_6_' + continent + '_ss_emt'
# Output directory
expodir = '/scratch/pj276/ifl_corridors/gfc/gfc1_6/gfc1_6_' + continent + '_ss_cs_' + scen
# Morphology directory
expmdir = '/scratch/pj276/ifl_corridors/gfc/gfc1_6/gfc1_6_morph'
# Morphology neighborhood (in # of cells)
mnbo = '33'
# Region - one of 'af', 'as', or 'sa'
reg = continent
# Shape parameter for exponential function.
# Use -0.01 to approximate a linear transformation from forest cover to suitability.
# Use -2 to approximate a non-linear transformation where suitability drops off steeply with forest cover decrease
# Use 2 to approximate a non-linear transformation where suitability drops off gradually with forest cover decrease
sdic = {sdic}
cshape = sdic[scen]
# Name corresponding to shape parameter
spname = scen
# Year
expy = expyear
# Use forest gain tree cover layer (only makes a difference if year is 2012)
usegain = 'no'
# Switch to use morphology or not
usemorph = 'no'

###############################################################################
# hf2_aggregate_rasters.py
###############################################################################
# Input dir
aggdir = '/scratch/pj276/ifl_corridors/gfc/gfc1_6/gfc1_6_' + continent + '_ss_cs_' + scen
# Output dir
aggodir = '/scratch/pj276/ifl_corridors/gfc/gfc1_6/gfc1_6_' + continent + '_ss_cs_' + scen + '_90m' 
# Amount by which to expand bounding box dimensions
exvec = [-90, -90, 90, 90]
# New cell dimension
aggdim = "90"
# No data value
ndval = "255"
# Source file name pattern
srcaggname = [continent + '_rtile_cs' + str(expy) + scen]
# Destination file name pattern
dstaggname = [continent + '_rtile_cs' + str(expy) + scen + '_90m']
"""

for i in ['af','as','sa']:
    for j in ['neg2', 'neg01', 'pos2']:
        for k in ['2000', '2012', '2018']:
            context = {"continent": i, "scen": j, "expyear": k, "sdic": {'neg2': -2.0, 'neg01': -0.01, 'pos2': 2.0}}
            with open(mdir + '/' + pbase + '_' + i + '_' + j + '_' + k + '.py','w') as myfile:
                myfile.write(ptwrite.format(**context))

#%%
#########################################################################
#---- hf2_make_exp_cs_submit_xx.sh ----
# Template for writing submit script
ttwrite = """#!/bin/bash
#SBATCH --job-name=mkcs_{continent}
#SBATCH --chdir=/scratch/pj276/ifl_corridors/output
#SBATCH --output=/scratch/pj276/ifl_corridors/logs/mkcs_%A_%a.log
#SBATCH --time=00:15:00
#SBATCH --partition=all
#SBATCH --cpus-per-task=1
#SBATCH --mem 24000
#SBATCH --array={tnums}

# 1-49 for the full set of tiles for south america
# 1-54 for the full set of tiles for africa
# 1-67 for the full set of tiles for asia

module load anaconda/latest
module load gdal
echo "File ID is:"

echo ${SLURM_ARRAY_TASK_ID}

srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepexp1_6/hf2_make_exp_cs_{continent}_{scen}_{expyear}.py ${SLURM_ARRAY_TASK_ID}
"""

# Get number of files in each continent folder to make submit scripts
# and write into the template and then write to file
for i in ['af','as','sa']:
    for j in ['neg2', 'neg01', 'pos2']:
        for k in ['2000', '2012', '2018']:
            # Edge matched tile directory
            d = gfcdir + '/' + os.path.basename(gfcdir) + '_' + i + '_ss_emt'
            # Get number of jobs for submit scripts using length of file list
            fl = os.listdir(d)
            fl = [z for z in fl if 'datamask' in z and z.endswith('.tif')]
            tnums = [np.int(z.split('_')[-1].split('.')[0]) for z in fl]
            tnums.sort()
            tnums = [str(z) for z in tnums]
            tnums = ["," . join(tnums)][0]
            # Set task ID variable (this is a nuisance, can't figure out how to ignore it in the template)
            SLURM_ARRAY_TASK_ID = 'SLURM_ARRAY_TASK_ID'
            # Dic that will be expanded in the format command
            context = {"tnums": tnums, "SLURM_ARRAY_TASK_ID": SLURM_ARRAY_TASK_ID,
                       "continent": i, "scen": j, "expyear": k}
            with open(mdir + '/hf2_make_exp_cs_submit_' + i + '_' + j + '_' + k + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))

#%%
#########################################################################
#---- Make hf2_make_exp_cs_xx.py scripts ----
for i in ['af','as','sa']:
    for j in ['neg2', 'neg01', 'pos2']:
        for k in ['2000', '2012', '2018']:
            copyfile(ssdir + '/hf2_make_exp_cs.py', mdir + '/hf2_make_exp_cs_'+ i + '_' + j + '_' + k + '.py')
            f = fileinput.FileInput(mdir + '/hf2_make_exp_cs_'+ i + '_' + j + '_' + k + '.py', inplace=True)
            for line in f:
                line = line.replace("import param_file as p1", "import " + pbase + '_' + i + '_' + j + '_' + k + " as p1")
                line = line.replace("sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')",
                                    "sys.path.append('/home/pj276/projects/ifl_corridors/code/scripts/prepexp1_6/')")
                print line,
            f.close()

##########################################################################
###---- SUBMIT JOBS ----
#for i in ['af', 'as', 'sa']:
#    for j in ['neg2', 'neg01', 'pos2']:
#        for k in ['2000', '2012', '2018']:
#            cmd = "sbatch " +   mdir + '/hf2_make_exp_cs_submit_' + i + '_' + j + '_' + k + '.sh'
#            print "Submitting Job " + i + '_' + j + '_' + k + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i + '_' + j + '_' + k

##########################################################################

#%%
#########################################################################
#---- hf2_aggregate_rasters_xx.sh ----
# Template for writing submit script
ttwrite = """#!/bin/bash
#SBATCH --job-name=aggrast
#SBATCH --chdir=/scratch/pj276/ifl_corridors/output
#SBATCH --output=/scratch/pj276/ifl_corridors/logs/aggrast_%A_%a.log
#SBATCH --time=00:60:00
#SBATCH --partition=all
#SBATCH --cpus-per-task=12
#SBATCH --mem 24000
module load anaconda/latest

srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepexp1_6/hf2_aggregate_rasters_{continent}_{scen}_{expyear}.py
"""
# Write template to file
for i in ['af','as','sa']:
    for j in ['neg2', 'neg01', 'pos2']:
        for k in ['2000', '2012', '2018']:
            # Dic that will be expanded in the format command
            context = {"continent": i, "scen": j, "expyear": k}
            with open(mdir + '/hf2_aggregate_rasters_submit_' + i + '_' + j + '_' + k + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))

#%%
#########################################################################
#---- Make hf2_aggregate_rasters_xx.py scripts ----
for i in ['af','as','sa']:
    for j in ['neg2', 'neg01', 'pos2']:
        for k in ['2000', '2012', '2018']:
            copyfile(ssdir + '/hf2_aggregate_rasters.py', mdir + '/hf2_aggregate_rasters_'+ i + '_' + j + '_' + k + '.py')
            f = fileinput.FileInput(mdir + '/hf2_aggregate_rasters_'+ i + '_' + j + '_' + k + '.py', inplace=True)
            for line in f:
                line = line.replace("import param_file as p1", "import " + pbase + '_' + i + '_' + j + '_' + k + " as p1")
                line = line.replace("sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')",
                                    "sys.path.append('/home/pj276/projects/ifl_corridors/code/scripts/prepexp1_6/')")
                print line,
            f.close()

##########################################################################
##---- SUBMIT JOBS ----
#for i in ['af', 'as', 'sa']:
#    for j in ['neg2', 'neg01', 'pos2']:
#        #for k in ['2000', '2012', '2018']:
#        for k in ['2018']:
#            cmd = "sbatch " +   mdir + '/hf2_aggregate_rasters_submit_' + i + '_' + j + '_' + k + '.sh'
#            print "Submitting Job " + i + '_' + j + '_' + k + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i + '_' + j + '_' + k
##########################################################################

