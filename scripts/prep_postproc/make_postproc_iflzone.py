# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 13:23:26 2020
Make scripts for making post processing scripts for ifl zones
@author: pj276
"""

#%% Modules
import commands, fileinput, os, subprocess, glob
from shutil import copyfile
import geopandas as gpd, rasterio
import numpy as np
import inspect
from shapely.geometry import box

#%%
# Monsoon git dir
gdir = '/home/pj276/projects/ifl_corridors'

# Source script directory
ssdir = gdir + '/code/scripts'

# Basename of parameter file
pbase = 'p1_012220_pp'

#%%
# Main directory holding param files and scripts
mdir = gdir + '/code/scripts/prep_postproc'
# Make dir if it doesn't exist already
if not os.path.exists(mdir):
    os.mkdir(mdir)

#%%
# Run switches
runP1 = True # make param files
runS1 = True # annual tree cover loss and flux calculations
runS2 = True # aggregate csvs produced by S1 or S2

#%%
# #########################################################################
if runP1 == True:
    #---- Section 1 parameter file ----
    # Template for parameter files (one for each region and exp scenario and year)
    ptwrite = """# -*- coding: utf-8 -*-
    \"\"\"Created on Sun Sep 29 20:54:49 2019
    This is a parameter file to make it easier to re-run corridor operations.
    @author: pj276
    \"\"\"
    
    # region code
    region = '{region}'
    # scenario
    scenario = '{scenario}'
    
    #----
    # Input dir
    idir = '/scratch/pj276/ifl_corridors'
    # Zone field
    field = 'gridcode'
    
    #----
    # Projection template
    rtemp = '/scratch/pj276/ifl_corridors/rp5k130k90m_' + scenario + '_' + region + '_2000_ifl/cfinal/flux/csum_rp5k130k90m_' + scenario + '_' + region + '_2000_ifl_mos.tif'
    # Input zone shapefile
    #inshp = '/scratch/pj276/undp_connectivity/vectors/WWF_ecoregions.shp'
    inshp = '/scratch/pj276/ifl_corridors/iflzones/{region}_iflzone_cs2000{scenario}_90m.shp'
    
    #----
    # Directory holding original flux files
    ofluxdir = 'flux'
    
    # Tree cover directory
    emtd = '/gfc/gfc1_6/gfc1_6_' + region + '_ss_emt_90m'
    
    # Output flux directory
    fluxdir = 'lybytc2000_iflzone'
    # use 'lybytc2000_iflzone' when calculating tree cover and flux
    # use 'lybytc2000_wwf' when calculating tree cover and flux
    # use 'flux_wwf' when calculating flux
    # use 'dflux_wwf' when calculating direct fluxes

    # Output flux name for aggregated csvs
    fluxcsv = 'iflzone_lybytc2000_flux_sums.csv'
    # use 'iflzone_lybytc2000_flux_sums.csv' when calculating tree cover and flux
    # use 'wwf_lybytc2000_flux_sums.csv' when calculating tree cover and flux
    # use 'wwf_flux_sums.csv' when calculating flux
    # use 'wwf_dflux_sums.csv' when calculating direct flux
    
    # Output flux directory for aggregated csvs
    aggcsvdir = 'iflzone_lybytc2000_fluxes'
    # use 'iflzone_lybytc2000_fluxes' when calculating tree cover and flux
    # use 'wwf_lybytc2000_fluxes' when calculating tree cover and flux
    # use 'wwf_fluxes' when calculating flux
    # use 'wwf_dfluxes' when calculating direct flux
    
    # Suffix to filter csv files
    endpatt = 'lybytc2000flux.csv'
    # use 'lybytc2000flux.csv' if calculating tree cover and flux
    # use 'flux.csv' if calculating fluxes
    # use 'fcincorrnoifl.csv' if calculating forest cover

    """
    # Fix extra tabs
    ptwrite = inspect.cleandoc(ptwrite)
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            context = {"region": i, "scenario": j, "sdic": {'neg2': -2.0, 'neg01': -0.01, 'pos2': 2.0}}
            with open(mdir + '/' + pbase + '_' + i + '_' + j + '.py', 'w') as myfile:
                myfile.write(ptwrite.format(**context))

#%%
def mainscriptmake(dictofreplacements, oldscriptdir, newscriptdir, **kwargs):
    """This function takes a dictionary of strings that are used to replace
    strings in an existing script. The new script is written to a new 
    directory. kwargs allows for iterated variables to be used as inputs.
    Generally, kwargs will be regions (e.g. Africa, Asia, South America),
    years (e.g. 2000, 2013), or scenarios (e.g. neg2, neg01, pos2) and are indexed
    as i, j, k, etc.
    """
    copyfile(oldscriptdir + '/' + dictofreplacements["scriptname"] + '.py', newscriptdir + '/' +
             dictofreplacements["scriptname"] + '_' + i + '_' + j + '.py')
    f = fileinput.FileInput(newscriptdir + '/' + dictofreplacements["scriptname"] + '_' + i + '_' + j + '.py',
                            inplace=True)
    for line in f:
        line = line.replace("import " + dictofreplacements["paramoldname"] + " as p1",
                            "import " + dictofreplacements["paramnewname"] + '_' + i + '_' + j + " as p1")
        line = line.replace("sys.path.append('" + dictofreplacements["syspathold"] + "')",
                            "sys.path.append('" + dictofreplacements["syspathnew"] + "')")
        print line,
    f.close()
    return(True)

#%%
##---- Make tile index needed by scripts
# Make tile index of flux files (use variables from above param file)
for i in ['af','as','sa']:
    for j in ['neg2', 'neg01', 'pos2']:
        for k in ['2000', '2013']:
            idir = '/scratch/pj276/ifl_corridors'
            atx = idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/flux/' + j + '_' + i + '_' + k + '_cflux.shp'
            if not os.path.exists(atx):
                path = idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/flux'
                df = gpd.GeoDataFrame(columns=['location','geometry'])
                for dir, subdir, files in os.walk(path):
                    for fname in files:
                        if fname.startswith('sumcorr') & fname.endswith(".tif"):
                            bounds = rasterio.open(os.path.join(dir +"/", fname)).bounds
                            df = df.append({'location':fname, 'geometry': box(bounds[0], bounds[1], bounds[2], bounds[3])},ignore_index=True)
                df.crs = rasterio.open(os.path.join(dir +"/", fname)).crs
                df.to_file(atx)

##---- Make vrts of annual tree cover area loss layers
for i in ['af','as','sa']:
    for k in range(0,19,1):
        tcvrt = idir + '/gfc/gfc1_6/gfc1_6_' + i + '_ss_emt_90m/' + i + '_' + 'lybytc2000_' + str(k) + '.vrt'
        if not os.path.exists(tcvrt):
            # List of files for vrt
            l1 = glob.glob(idir + '/gfc/gfc1_6/gfc1_6_' + i + '_ss_emt_90m/*')
            l1 = [p for p in l1 if 'lybytc2000_' + str(k) + '_' in p and p.endswith('.tif') and 'mos' not in p]
            # Build a vrt from a file list
            acmd = ["gdalbuildvrt", tcvrt]
            for myfile in l1:
                acmd.append(myfile)
            subprocess.call(acmd)

##---- Make vrts of tree cover 2000 layers
for i in ['af','as','sa']:
    tcvrt = idir + '/gfc/gfc1_6/gfc1_6_' + i + '_ss_emt_90m/' + i + '_' + 'treecover2000.vrt'
    if not os.path.exists(tcvrt):
        # List of files for vrt
        l1 = glob.glob(idir + '/gfc/gfc1_6/gfc1_6_' + i + '_ss_emt_90m/*')
        l1 = [p for p in l1 if 'treecover2000' in p and p.endswith('.tif') and 'mos' not in p]
        # Build a vrt from a file list
        acmd = ["gdalbuildvrt", tcvrt]
        for myfile in l1:
            acmd.append(myfile)
        subprocess.call(acmd)

##---- Make directories needed by scripts if they don't exist
for i in ['af','as','sa']:
    for j in ['neg2', 'neg01', 'pos2']:
        for k in ['2000', '2013']:
            idir = '/scratch/pj276/ifl_corridors'
            # Copy directory name from fluxdir above
            if not os.path.exists(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/lybytc2000_iflzone'):
                os.mkdir(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/lybytc2000_iflzone')

# Get directory names from aggcsvdir above
if runS1 == True and runS2 == True:
    if not os.path.exists(idir + '/iflzone_lybytc2000_fluxes'):
        os.mkdir(idir + '/iflzone_lybytc2000_fluxes')

#%%
if runS1 == True:
    ##---- Make hf2_lyby2000_iflzone_stats_xx.py scripts ----
    repldict = {"scriptname": "hf2_lyby2000_iflzone_stats", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prep_postproc/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_lyby2000_iflzone_stats_submit_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=lyby2ks
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/lyby2ks_%A_%a.log
    #SBATCH --time=00:03:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=2
    #SBATCH --mem 4000
    #SBATCH --array={njsfupdate}-{njefupdate}
    
    # Number of tasks should range from 1 to the number of shapes in 
    # the zone.
    
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prep_postproc/hf2_lyby2000_iflzone_stats_{region}_{scenario}.py ${SLURM_ARRAY_TASK_ID}
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Get number of jobarrays
    # and write into the template and then write to file
    
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            inshp = '/scratch/pj276/ifl_corridors/iflzones/' + i + '_iflzone_cs2000' + j + '_90m.shp'
            if os.path.exists(inshp):
                shpfile = gpd.read_file(inshp)
                ucode = np.unique(shpfile['gridcode'].tolist())
                myjobcount =  len(ucode)

                # Get number of jobs
                njsfupdate = 1
                njefupdate = myjobcount
                # Set task ID variable (this is a nuisance, can't figure out how to ignore it in the template)
                SLURM_ARRAY_TASK_ID = 'SLURM_ARRAY_TASK_ID'
                
                context = {
                 "njsfupdate": njsfupdate, 
                 "njefupdate": njefupdate,
                 "SLURM_ARRAY_TASK_ID": SLURM_ARRAY_TASK_ID,
                 "region": i,
                 "scenario": j
                }
                with open(mdir + '/hf2_lyby2000_iflzone_stats_submit_' + i + '_' + j + '.sh','w') as myfile:
                    myfile.write(ttwrite.format(**context))
        
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            cmd = "sbatch " +   mdir + '/hf2_lyby2000_iflzone_stats_submit_' + i + '_' + j + '.sh'
#            print "Submitting Job " + i + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i

#%%
if runS2 == True:
    ##---- Make hf2_aggregate_flux_csvs_xx.py scripts ----
    repldict = {"scriptname": "hf2_aggregate_flux_csvs", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prep_postproc/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
            
    #---- hf2_aggregate_flux_csvs_submit_xx.sh ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=prpcsvs
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/prpcsvs_%A_%a.log
    #SBATCH --time=00:03:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=1
    #SBATCH --mem 4000
    
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prep_postproc/hf2_aggregate_flux_csvs_{region}_{scenario}.py
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Write submit script to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            context = {"region": i, "scenario": j}
            with open(mdir + '/hf2_aggregate_flux_csvs_submit_' + i + '_' + j + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))        
    
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            cmd = "sbatch " +   mdir + '/hf2_aggregate_flux_csvs_submit_' + i + '_' + j + '.sh'
#            print "Submitting Job " + i + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i