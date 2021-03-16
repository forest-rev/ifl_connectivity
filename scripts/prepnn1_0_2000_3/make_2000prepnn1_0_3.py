# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 14:51:46 2020
Make scripts for making point pair scripts
@author: pj276
"""

#%% Modules
import commands, fileinput, os
from shutil import copyfile
import inspect

#%%
# Monsoon git dir
gdir = '/home/pj276/projects/ifl_corridors'

# Source script directory
ssdir = gdir + '/code/scripts'

# Basename of parameter file
pbase = 'p1_030120_pp'

#%%
# Main directory holding param files and scripts
mdir = gdir + '/code/scripts/prepnn1_0_2000_3'
# Make dir if it doesn't exist already
if not os.path.exists(mdir):
    os.mkdir(mdir)

#%%
# Run switches
runP1 = True # make param files

runS1 = False # hf2_make_directories_xx.py
runS2 = False # hf2_point_pair_create
runS3 = False # hf2_point_pair_create_concat
runS4 = False # hf2_get_line_count_for_ppc2
runS5 = False # hf2_point_pair_create_part2 Check file counts.
runS6 = False # hf2_point_pair_create_part3  
runS7 = False # hf2_point_pair_create_part4
runS8 = False # hf2_ppc4_redo_list This may not work properly. Check file counts.
# May have to rerun S7 at this point
runS9 = False # hf2_ptnbs_rename
runS10 = False # hf2_point_pair_create_part5
runS11 = False # hf2_get_line_count_for_correx
runS12 = False # hf2_pp_corridor_expander_alt3
runS13 = False # hf2_get_corridor_unique_id_count
runS14 = True # hf2_sum_level1_corridors
runS15 = False # hf2_clean_up
runS16 = False # hf2_summarize_corridors_alt
runS17 = False # hf2_stack_corridors
runS18 = False # hf2_mosaic_corridor_tiles
runS19 = False # hf2_aggregate_rasters
runS20 = False # hf2_make_convex_hulls
runS21 = False # hf2_get_convex_hull_areas

#%%
# #########################################################################
if runP1 == True:
    #---- Section 1 parameter file ----
    # Template for parameter files (one for each region and exp scenario and year)
    ptwrite = """# -*- coding: utf-8 -*-
    \"\"\"Created on Thu Jun 4 10:20:00 2019
    This is a parameter file to make it easier to re-run corridor operations.
    @author: pj276
    \"\"\"
    
    # region code
    region = '{region}'
    # scenario
    scenario = '{scenario}'
    
    #----
    # Check wdir, inpoly, fnet, and year
    wdir = 'rp5k130k90m_' + scenario + '_' + region + '_2000_ifl'
    fpre = 'ptnbs_130k'
    #----
    
    #----
    # hf2_make_directories.py
    dlist = ['/distfiles/ptile', '/cfinal/flux', '/cfinal/nprob', '/cfinal/dstatsl1', '/cfinal/dstatsl2', '/clevel1count', '/corridors', '/finfiles', '/ppc4redo', '/convex_hulls', '/pointfiles', '/areadfs']
    #----
    
    #----
    # hf2_point_pair_create.py
    # Set up inputs
    odir = '/scratch/pj276/ifl_corridors/' + wdir
    # Point directory
    pointdir = '/pointfiles'
    xspac = 5000 # x spacing in meters
    yspac = 5000 # y spacing in meters
    
    # Cell size (interval)
    xint = 90 # cell size
    yint = 90 # cell size
    
    #ss
    if region == 'af':
        # af
        snap_x = -5009557.085699999704957 + (10000*xint)
        snap_y = 2212524.833229999989271 - (10000*yint)
    if region == 'as':
        #as
        snap_x = -5009551.383070000447333 + (10000*xint)
        snap_y = 3320266.254170000087470 - (10000*yint)
    if region == 'sa':
        #sa
        snap_x = -5009557.085699999704957 + (10000*xint)
        snap_y = 2212524.833229999989271 - (10000*yint)
    

    fieldname = 'GRIDCODE'
    tdist = 130000 # maximum point distance to consider, in meters
    dname = '5k' # code to identify point type (regular point)
    ##inpoly = '/projects/above_gedi/pjantz/ifl_corridors_data/hifo/hland_tropics_ssu_' + region + '.shp'
    inpoly = '/scratch/pj276/ifl_corridors/ifl/ifl_2000_tropics_ssu_' + region + '.shp'
    csvout =  '/scratch/pj276/ifl_corridors/' + wdir + '/' + fpre + '.csv'
    pointlog = '/scratch/pj276/ifl_corridors/' + wdir + '/' + 'pointlog.csv'
    #----
    
    #----
    # hf2_point_pair_create_concat.py
    # Nothing needed
    #----
    
    #----
    #hf2_get_line_count_for_ppc2.py
    outcount = '/scratch/pj276/ifl_corridors/' + wdir + '/ppc2count.txt'
    
    #----
    # hf2_point_pair_create_part2.py
    nbcsvpre = '/scratch/pj276/ifl_corridors/' + wdir + '/' + fpre + '_'
    #----
    
    #----
    # hf2_point_pair_create_part3.py
    nbcsvall = '/scratch/pj276/ifl_corridors/' + wdir + '/' + fpre + '_all.csv'                                                                                                             
    spat = '/scratch/pj276/ifl_corridors/' + wdir + '/' + fpre + '_*'
    ppcount = '/scratch/pj276/ifl_corridors/' + wdir + '/ppcount'
    rdfsub = '/scratch/pj276/ifl_corridors/' + wdir + '/rdfsub'
    #----
    
    #----
    # hf2_point_pair_create_part4.py
    # Directory to hold redo files (directory is created in part3)
    ppc4redodir = '/ppc4redo'
    # Convex hull directory
    chulldir = '/convex_hulls'
    # Directory holding area tables
    areadfs = '/areadfs'
    # Area data frame name
    adfname_temp = odir + areadfs + '/areadf_'
    # Tree cover vrt (year should correspond with cost surface)
    tcvrt = '/scratch/pj276/ifl_corridors/gfc/gfc1_6/gfc1_6_' + region + '_ss_emt_90m/mos_tc2000.vrt'
    
    #----
    
    #----
    # hf2_ppc4_redo_list.py
    ppc4jnums = None
    ppc4rdtxt = odir + '/redolist.txt'
    
    #----
    
    #----
    # hf2_ptnbs_rename.py
    srcname1 = odir + '/ptnbs_130k_all.csv'
    dstname1 = odir + '/ptnbs_130k_all_orig.csv'
    srcname2 = odir + '/rdfsub.csv'
    dstname2 = odir + '/ptnbs_130k_all.csv'
    
    #----
    
    #----
    # hf2_point_pair_create_part5.py
    adfname = odir + '/areadf.csv'
    
    #----
    
    #----
    # hf2_get_line_count_for_correx.py
    outcexcount = odir + '/correxcount.txt'
    emd = '/gfc/gfc1_6/gfc1_6_' + region + '_ss_cs_' + scenario + '_90m'
    forfileindex = '/scratch/pj276/ifl_corridors' + emd + '/cs2000' + scenario + 'index.txt'    
    
    #----
    
    #----
    # hf2_pp_corridor_expander_alt3.py
    ffdir = '/finfiles'
    sdir = '/scratch/pj276/ifl_corridors/' + wdir + '/corridors'
    # Input csv folder and base name
    icsv = '/scratch/pj276/ifl_corridors/' + wdir + '/' + fpre + '_all.csv'
    # File name for unfinished corridors for re-running
    # Created in another error logging script
    icsv2 = '/scratch/pj276/ifl_corridors/' + wdir + '/' + fpre + '_all_redfile.csv'
    # Name of redo file
    rdname = 'redome.txt'
    # vrt file prefix
    vrtp = '/cs2000_' + scenario + '_'
    # output corridor prefix
    ocp = '/corr_' + scenario + '_'
    # Output stats path names
    olcpp = '/distfiles/lcp'
    optilep = '/distfiles/ptile'
    xthp = 10 # percentile threshold to use for delineating corridors

    #----
       
    #----
    # hf2_get_corridor_unique_id_count.py
    finstatusptn = '/finstatus_*'
    outidcount = '/scratch/pj276/ifl_corridors/' + wdir + '/pidcount.txt'
    findex = '/scratch/pj276/ifl_corridors/' + wdir + '/findex.csv'
    ucpfile = '/scratch/pj276/ifl_corridors/' + wdir + '/ucpid.csv'
    
    #----
        
    #----
    # hf2_sum_level1_corridors.py
    ayear = '2000'
    ccdirp = '/clevel1count'

    #----
    
    #----
    # hf2_summarize_corridors_alt.py
    # Csv file prefix
    csvfp = 'lcc_ptile_dists_' + ayear + '_' + scenario + '_'
    # Output stats csv basename for all stat values per patch pair
    statsoutl1 = odir + '/cfinal/dstatsl1/dist_area_statsl1_'
    # Output stats csv basename for single stat value per patch pair
    statsoutl2 = odir + '/cfinal/dstatsl2/dist_area_statsl2_'
    # Get spatial reference
    prjfile = '/scratch/pj276/ifl_corridors/prjfiles/ss_' + region + '.prj'
    # Get interaction distance (median dispersal distance)
    ixdist = 25

    #----
    
    #----
    # hf2_stack_corridors.py
    fnet = '/scratch/pj276/ifl_corridors/fishnets/fnet_' + region + '_270k_ss.shp'

    #----
    
    #----
    # hf2_aggregate_rasters.py
    aggdir = odir + '/cfinal/flux'    
    # Output dir
    aggodir = odir + '/cfinal/flux'    
    # The amount by which to expand bounding box dimensions
    exvec = [-540, -540, 540, 540]    
    # Cell dim
    aggdim = "540"    
    # No data val
    ndval = "255"        
    # Source file name pattern
    srcaggname = ['csum_rp5k130k90m_' + scenario + '_' + region + '_' + ayear + '_ifl_mos']
    # Replacement file name pattern
    dstaggname = ['csum_rp5k130k540m_' + scenario + '_' + region + '_' + ayear + '_ifl_mos']
    
    #----

    """
    # Fix extra tabs
    ptwrite = inspect.cleandoc(ptwrite)
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            context = {"region": i, "scenario": j, "sdic": {'neg2': -2.0, 'neg01': -0.01, 'pos2': 2.0}}
            with open(mdir + '/' + pbase + '_' + i + '_' + j + '.py', 'w') as myfile:
                myfile.write(ptwrite.format(**context))

#%%
#########################################################################
# MAIN SCRIPTS
#########################################################################
            
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
if runS1 == True:
    ##---- Make hf2_make_directories_xx.py scripts ----
    repldict = {"scriptname": "hf2_make_directories", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
            
    #---- hf2_make_directories_submit_xx.sh ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=mkdir
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/mkdir_%A_%a.log
    #SBATCH --time=00:00:60
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=1
    #SBATCH --mem 4000
    
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_make_directories_{region}_{scenario}.py
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Write submit script to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            context = {"region": i, "scenario": j}
            with open(mdir + '/hf2_make_directories_submit_' + i + '_' + j + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))        
    
    #---- SUBMIT JOBS ----
#    S1jnumlist = []
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            cmd = "sbatch " +   mdir + '/hf2_make_directories_submit_' + i + '_' + j + '.sh'
#            print "Submitting Job " + i + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            S1jnumlist.append(jobnum.split(' ')[-1])
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i + '_' + j

#%%
if runS2 == True:
    ##---- Make hf2_point_pair_create_xx.py scripts ----
    repldict = {"scriptname": "hf2_point_pair_create", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_point_pair_create_submit_xx.sh scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=ppcmap
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/ppcmap_%A_%a.log
    #SBATCH --time=00:30:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=1
    #SBATCH --mem 16000
    
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_point_pair_create_{region}_{scenario}.py
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Write submit script to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            context = {"region": i, "scenario": j}
            with open(mdir + '/hf2_point_pair_create_submit_' + i + '_' + j + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))        
    
    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            cmd = "sbatch --depend=afterok:{0}:{1} ".format(S1jnumlist[0],S1jnumlist[-1]) + mdir + '/hf2_point_pair_create_submit_' + i + '_' + j + '.sh'
#            print "Submitting Job " + i + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i

#%%
if runS3 == True:
    ##---- Make hf2_point_pair_create_concat_xx.py scripts ----
    repldict = {"scriptname": "hf2_point_pair_create_concat", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_point_pair_create_concat_submit_xx.sh scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=ppcat
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/ppcat_%A_%a.log
    #SBATCH --time=00:03:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=1
    #SBATCH --mem 4000
    
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_point_pair_create_concat_{region}_{scenario}.py
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Write submit script to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            context = {"region": i, "scenario": j}
            with open(mdir + '/hf2_point_pair_create_concat_submit_' + i + '_' + j + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))      
    
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            cmd = "sbatch " +   mdir + '/hf2_point_pair_create_concat_submit_' + i + '_' + j + '.sh'
#            print "Submitting Job " + i + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i

#%%
if runS4 == True:
    ##---- Make hf2_get_line_count_for_ppc2_xx.py scripts ----
    repldict = {"scriptname": "hf2_get_line_count_for_ppc2", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_get_line_count_for_ppc2_submit_xx.sh scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=ppc2count
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/ppc2c_%A_%a.log
    #SBATCH --time=00:00:30
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=1
    #SBATCH --mem 4000
    
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_get_line_count_for_ppc2_{region}_{scenario}.py
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Write submit script to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            context = {"region": i, "scenario": j}
            with open(mdir + '/hf2_get_line_count_for_ppc2_submit_' + i + '_' + j + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))
    
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            cmd = "sbatch " +   mdir + '/hf2_get_line_count_for_ppc2_submit_' + i + '_' + j + '.sh'
#            print "Submitting Job " + i + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i

#%%
if runS5 == True:
    ##---- Make hf2_point_pair_create_part2_xx.py scripts ----
    repldict = {"scriptname": "hf2_point_pair_create_part2", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_point_pair_create_part2_submit_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=ppcmap2
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/ppcm_%A_%a.log
    #SBATCH --time=00:03:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=1
    #SBATCH --mem 16000
    #SBATCH --array={njsfupdate}-{njefupdate}
    
    # The task range should range from 1 to the number of entries in the ptnbs_xx_xx.csv
    # divided by 10 with a ceiling operator. Should be given by ppc2count.txt.
    
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_point_pair_create_part2_{region}_{scenario}.py ${SLURM_ARRAY_TASK_ID}
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Get number of jobarrays
    # and write into the template and then write to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
            odir = '/scratch/pj276/ifl_corridors/' + wdir
            if os.path.exists(odir + '/ppc2count.txt'):
                with open(odir + '/ppc2count.txt', 'r') as countfile:
                    myjobcount = int(float(countfile.readline()))
                # Get number of jobs for update forest cover script using length of file list
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
                with open(mdir + '/hf2_point_pair_create_part2_submit_' + i + '_' + j + '.sh','w') as myfile:
                    myfile.write(ttwrite.format(**context))
    
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
#            odir = '/scratch/pj276/ifl_corridors/' + wdir
#            if os.path.exists(odir + '/ppc2count.txt'):
#                cmd = "sbatch " +   mdir + '/hf2_point_pair_create_part2_submit_' + i + '_' + j + '.sh'
#                print "Submitting Job " + i + " with command: %s" % cmd
#                status, jobnum = commands.getstatusoutput(cmd)
#                if (status == 0):
#                    print "Job1 is %s" % jobnum
#                else:
#                    print "Error submitting Job " + i

#%%
if runS6 == True:
    ##---- Make hf2_point_pair_create_part3_xx.py scripts ----
    repldict = {"scriptname": "hf2_point_pair_create_part3", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_point_pair_create_part3_submit_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=ppcp3
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/ppcp3_%A_%a.log
    #SBATCH --time=00:25:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=1
    #SBATCH --mem 32000
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_point_pair_create_part3_{region}_{scenario}.py
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Write submit script to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            context = {"region": i, "scenario": j}
            with open(mdir + '/hf2_point_pair_create_part3_submit_' + i + '_' + j + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))
    
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            cmd = "sbatch " +   mdir + '/hf2_point_pair_create_part3_submit_' + i + '_' + j + '.sh'
#            print "Submitting Job " + i + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i

#%%
if runS7 == True:
    ##---- Make hf2_point_pair_create_part4_xx.py scripts ----
    repldict = {"scriptname": "hf2_point_pair_create_part4", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_point_pair_create_part4_submit_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=ppcp4
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/ppcp4_%A_%a.log
    #SBATCH --time=00:30:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=2
    #SBATCH --mem 16000
    #SBATCH --array={njsfupdate}-{njefupdate}
    
    # I think the number of jobs is the number of rows in ppcount.csv minus 1.
    # For big jobs, use 1.5 hours
    
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_point_pair_create_part4_{region}_{scenario}.py ${SLURM_ARRAY_TASK_ID}
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Get number of jobarrays
    # and write into the template and then write to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
            odir = '/scratch/pj276/ifl_corridors/' + wdir
            if os.path.exists(odir + '/ppcount.csv'):
                with open(odir + '/ppcount.csv', 'r') as countfile:
                    myjobcount = sum(1 for line in countfile)-1
                # Get number of jobs for update forest cover script using length of file list
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
                with open(mdir + '/hf2_point_pair_create_part4_submit_' + i + '_' + j + '.sh','w') as myfile:
                    myfile.write(ttwrite.format(**context))
    
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
#            odir = '/scratch/pj276/ifl_corridors/' + wdir
#            if os.path.exists(odir + '/ppcount.csv'):
#                cmd = "sbatch " +   mdir + '/hf2_point_pair_create_part4_submit_' + i + '_' + j + '.sh'
#                print "Submitting Job " + i + " with command: %s" % cmd
#                status, jobnum = commands.getstatusoutput(cmd)
#                if (status == 0):
#                    print "Job1 is %s" % jobnum
#                else:
#                    print "Error submitting Job " + i

#%%
if runS8 == True:
    ##---- Make hf2_ppc4_redo_list_xx.py scripts ----
    repldict = {"scriptname": "hf2_ppc4_redo_list", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Update parameter file with job count
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
            odir = '/scratch/pj276/ifl_corridors/' + wdir
            if os.path.exists(odir + '/ppcount.csv'):
                with open(odir + '/ppcount.csv', 'r') as countfile:
                    myjobcount = sum(1 for line in countfile)-1
                f = fileinput.FileInput(mdir + '/' + pbase + '_' + i + '_' + j + '.py', inplace=True)
                for line in f:
                    line = line.replace("ppc4jnums = None", "ppc4jnums = " + str(myjobcount))
                    print line,
    
    ##---- Make hf2_ppc4_redo_list_submit_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=ppc4r
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/ppc4r_%A_%a.log
    #SBATCH --time=00:03:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=1
    #SBATCH --mem 4000
    
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_ppc4_redo_list_{region}_{scenario}.py
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Write submit script to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            context = {"region": i, "scenario": j}
            with open(mdir + '/hf2_ppc4_redo_list_submit_' + i + '_' + j + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))
    
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            cmd = "sbatch " +   mdir + '/hf2_ppc4_redo_list_submit_' + i + '_' + j + '.sh'
#            print "Submitting Job " + i + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i
#    
#    # Check if there are files to rerun
#    import sys, time
#    rerunflag = []
#    for i in ['af','as','sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
#            odir = '/scratch/pj276/ifl_corridors/' + wdir
#            while not os.path.exists(odir + '/redolist.txt'):
#               time.sleep(5) 
#            with open(odir + '/redolist.txt', 'r') as rdfile:
#                mylines = [myline for myline in rdfile]
#            if len(mylines) > 0 and mylines[0] == 'None':
#                print 'No files to rerun'
#            elif len(mylines) > 0:
#                print 'Rerun point pair create part 4 for ' + i + ' ' + j
#                rerunflag.append(i + ' ' + j)
#    if len(rerunflag) > 0:
#        sys.exit("Rerun point pair create part 4")

#%%
if runS9 == True:
    ##---- Make hf2_ptnbs_rename_xx.py scripts ----
    repldict = {"scriptname": "hf2_ptnbs_rename", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_ptnbs_rename_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=ptrname
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/ptrname_%A_%a.log
    #SBATCH --time=00:01:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=1
    #SBATCH --mem 4000
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_ptnbs_rename_{region}_{scenario}.py
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Write submit script to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            context = {"region": i, "scenario": j}
            with open(mdir + '/hf2_ptnbs_rename_submit_' + i + '_' + j + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))
    
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            cmd = "sbatch " +   mdir + '/hf2_ptnbs_rename_submit_' + i + '_' + j + '.sh'
#            print "Submitting Job " + i + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i

#%%
if runS10 == True:
    ##---- Make hf2_point_pair_create_part5_xx.py scripts ----
    repldict = {"scriptname": "hf2_point_pair_create_part5", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_point_pair_create_part5_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=ppcp5
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/ppcp5_%A_%a.log
    #SBATCH --time=00:20:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=1
    #SBATCH --mem 16000
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_point_pair_create_part5_{region}_{scenario}.py
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Write submit script to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            context = {"region": i, "scenario": j}
            with open(mdir + '/hf2_point_pair_create_part5_submit_' + i + '_' + j + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))
    
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            cmd = "sbatch " +   mdir + '/hf2_point_pair_create_part5_submit_' + i + '_' + j + '.sh'
#            print "Submitting Job " + i + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i

#%%
if runS11 == True:
    ##---- Make hf2_get_line_count_for_correx_xx.py scripts ----
    repldict = {"scriptname": "hf2_get_line_count_for_correx", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_get_line_count_for_correx_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=glcc
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/glcc_%A_%a.log
    #SBATCH --time=00:05:30
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=1
    #SBATCH --mem 8000
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_get_line_count_for_correx_{region}_{scenario}.py
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Write submit script to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            context = {"region": i, "scenario": j}
            with open(mdir + '/hf2_get_line_count_for_correx_submit_' + i + '_' + j + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))
    
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            cmd = "sbatch " +   mdir + '/hf2_get_line_count_for_correx_submit_' + i + '_' + j + '.sh'
#            print "Submitting Job " + i + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i
            

#%%
if runS12 == True:
    ##---- Make hf2_pp_corridor_expander_alt3_xx.py scripts ----
    repldict = {"scriptname": "hf2_pp_corridor_expander_alt3", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_pp_corridor_expander_alt3_submit_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=ppcex
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/ppcex_%A_%a.log
    #SBATCH --time=01:15:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=2
    #SBATCH --mem 4000
    #SBATCH --array={njsfupdate}-{njefupdate}

    module load anaconda/latest
    module load gdal
    echo "File ID is:"
    
    echo ${SLURM_ARRAY_TASK_ID}
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_pp_corridor_expander_alt3_{region}_{scenario}.py ${SLURM_ARRAY_TASK_ID}
      
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Get number of jobarrays
    # and write into the template and then write to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
            odir = '/scratch/pj276/ifl_corridors/' + wdir
            if os.path.exists(odir + '/correxcount.txt'):
                with open(odir + '/correxcount.txt', 'r') as countfile:
                    myjobcount = int(float(countfile.readline()))
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
                with open(mdir + '/hf2_pp_corridor_expander_alt3_submit_' + i + '_' + j + '.sh','w') as myfile:
                    myfile.write(ttwrite.format(**context))
    
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
#            odir = '/scratch/pj276/ifl_corridors/' + wdir
#            if os.path.exists(odir + '/correxcount.txt'):
#                cmd = "sbatch " +   mdir + '/hf2_pp_corridor_expander_alt3_submit_' + i + '_' + j + '.sh'
#                print "Submitting Job " + i + " with command: %s" % cmd
#                status, jobnum = commands.getstatusoutput(cmd)
#                if (status == 0):
#                    print "Job1 is %s" % jobnum
#                else:
#                    print "Error submitting Job " + i

#%%
if runS13 == True:
    ##---- Make hf2_get_corridor_unique_id_count_xx.py scripts ----
    repldict = {"scriptname": "hf2_get_corridor_unique_id_count", "paramoldname": "param_file", "paramnewname": pbase,
                "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
                "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_get_corridor_unique_id_count_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=cuid
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/cuid_%A_%a.log
    #SBATCH --time=00:30:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=1
    #SBATCH --mem 4000
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_get_corridor_unique_id_count_{region}_{scenario}.py
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Write submit script to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            context = {"region": i, "scenario": j}
            with open(mdir + '/hf2_get_corridor_unique_id_count_submit_' + i + '_' + j + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))

#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            cmd = "sbatch " +   mdir + '/hf2_get_corridor_unique_id_count_submit_' + i + '_' + j + '.sh'
#            print "Submitting Job " + i + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i

#%%
if runS14 == True:
    ##---- Make hf2_sum_level1_corridors_xx.py scripts ----
    repldict = {"scriptname": "hf2_sum_level1_corridors", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_sum_level1_corridors_submit_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=sl1
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/sl1_%A_%a.log
    #SBATCH --time=00:02:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=2
    #SBATCH --mem 8000
    #SBATCH --array={njsfupdate}-{njefupdate}

    module load anaconda/latest

    echo "File ID is:"
    
    echo ${SLURM_ARRAY_TASK_ID}
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_sum_level1_corridors_{region}_{scenario}.py ${SLURM_ARRAY_TASK_ID}
      
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Get number of jobarrays
    # and write into the template and then write to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
            odir = '/scratch/pj276/ifl_corridors/' + wdir
            if os.path.exists(odir + '/pidcount.txt'):
                with open(odir + '/pidcount.txt', 'r') as countfile:
                    myjobcount = int(float(countfile.readline()))
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
                with open(mdir + '/hf2_sum_level1_corridors_submit_' + i + '_' + j + '.sh','w') as myfile:
                    myfile.write(ttwrite.format(**context))
    
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
#            odir = '/scratch/pj276/ifl_corridors/' + wdir
#            if os.path.exists(odir + '/pidcount.txt'):
#                cmd = "sbatch " +   mdir + '/hf2_sum_level1_corridors_submit_' + i + '_' + j + '.sh'
#                print "Submitting Job " + i + " with command: %s" % cmd
#                status, jobnum = commands.getstatusoutput(cmd)
#                if (status == 0):
#                    print "Job1 is %s" % jobnum
#                else:
#                    print "Error submitting Job " + i

#%%
if runS15 == True:
    ##---- Make hf2_clean_up_xx.py scripts ----
    repldict = {"scriptname": "hf2_clean_up", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_clean_up_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=clup
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/clup_%A_%a.log
    #SBATCH --time=00:06:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=1
    #SBATCH --mem 8000
    
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_clean_up_{region}_{scenario}.py
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Write submit script to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            context = {"region": i, "scenario": j}
            with open(mdir + '/hf2_clean_up_submit_' + i + '_' + j + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))

#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            cmd = "sbatch " +   mdir + '/hf2_clean_up_submit_' + i + '_' + j + '.sh'
#            print "Submitting Job " + i + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i

#%%
if runS16 == True:
    ##---- Make hf2_summarize_corridors_alt_xx.py scripts ----
    repldict = {"scriptname": "hf2_summarize_corridors_alt", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_summarize_corridor_alt_submit_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=scorr
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/scorr_%A_%a.log
    #SBATCH --time=00:02:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=2
    #SBATCH --mem 8000
    #SBATCH --array={njsfupdate}-{njefupdate}

    module load anaconda/latest

    echo "File ID is:"
    
    echo ${SLURM_ARRAY_TASK_ID}
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_summarize_corridors_alt_{region}_{scenario}.py ${SLURM_ARRAY_TASK_ID}
      
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Get number of jobarrays
    # and write into the template and then write to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
            odir = '/scratch/pj276/ifl_corridors/' + wdir
            if os.path.exists(odir + '/pidcount.txt'):
                with open(odir + '/pidcount.txt', 'r') as countfile:
                    myjobcount = int(float(countfile.readline()))
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
                with open(mdir + '/hf2_summarize_corridors_alt_submit_' + i + '_' + j + '.sh','w') as myfile:
                    myfile.write(ttwrite.format(**context))
    
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
#            odir = '/scratch/pj276/ifl_corridors/' + wdir
#            if os.path.exists(odir + '/pidcount.txt'):
#                cmd = "sbatch " +   mdir + '/hf2_summarize_corridors_alt_submit_' + i + '_' + j + '.sh'
#                print "Submitting Job " + i + " with command: %s" % cmd
#                status, jobnum = commands.getstatusoutput(cmd)
#                if (status == 0):
#                    print "Job1 is %s" % jobnum
#                else:
#                    print "Error submitting Job " + i

#%%
if runS17 == True:
    ##---- Make hf2_stack_corridors_xx.py scripts ----
    repldict = {"scriptname": "hf2_stack_corridors", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_stack_corridors_submit_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=cstack
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/cstack_%A_%a.log
    #SBATCH --time=00:20:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=2
    #SBATCH --mem 64000
    #SBATCH --array={njsfupdate}-{njefupdate}

    module load anaconda/latest

    echo "File ID is:"
    
    echo ${SLURM_ARRAY_TASK_ID}
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_stack_corridors_{region}_{scenario}.py ${SLURM_ARRAY_TASK_ID}
      
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Get number of jobarrays
    # and write into the template and then write to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
            odir = '/scratch/pj276/ifl_corridors/' + wdir
            if i == 'af':
                # Set number of jobs 
                njsfupdate = 1
                njefupdate = 693
            elif i == 'as':
                # Set number of jobs 
                njsfupdate = 1
                njefupdate = 1100
            else:
                # Set number of jobs 
                njsfupdate = 1
                njefupdate = 676
            # Set task ID variable (this is a nuisance, can't figure out how to ignore it in the template)
            SLURM_ARRAY_TASK_ID = 'SLURM_ARRAY_TASK_ID'
            context = {
             "njsfupdate": njsfupdate, 
             "njefupdate": njefupdate,
             "SLURM_ARRAY_TASK_ID": SLURM_ARRAY_TASK_ID,
             "region": i,
             "scenario": j
            }
            with open(mdir + '/hf2_stack_corridors_submit_' + i + '_' + j + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))
    
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
#            odir = '/scratch/pj276/ifl_corridors/' + wdir            
#            cmd = "sbatch " +   mdir + '/hf2_stack_corridors_submit_' + i + '_' + j + '.sh'
#            print "Submitting Job " + i + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i

#%%
if runS18 == True:
    ##---- Make hf2_mosaic_corridor_tiles_xx.py scripts ----
    repldict = {"scriptname": "hf2_mosaic_corridor_tiles", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_mosaic_corridor_tiles_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=cmos
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/cmos_%A_%a.log
    #SBATCH --time=00:08:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=1
    #SBATCH --mem 16000
    
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_mosaic_corridor_tiles_{region}_{scenario}.py
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Write submit script to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            context = {"region": i, "scenario": j}
            with open(mdir + '/hf2_mosaic_corridor_tiles_submit_' + i + '_' + j + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))

#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            cmd = "sbatch " +   mdir + '/hf2_mosaic_corridor_tiles_submit_' + i + '_' + j + '.sh'
#            print "Submitting Job " + i + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i

#%%
if runS19 == True:
    ##---- Make hf2_aggregate_rasters_xx.py scripts ----
    repldict = {"scriptname": "hf2_aggregate_rasters", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_aggregate_rasters_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=ragg
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/ragg_%A_%a.log
    #SBATCH --time=00:10:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=12
    #SBATCH --mem 16000
    
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_aggregate_rasters_{region}_{scenario}.py
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Write submit script to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            context = {"region": i, "scenario": j}
            with open(mdir + '/hf2_aggregate_rasters_submit_' + i + '_' + j + '.sh','w') as myfile:
                myfile.write(ttwrite.format(**context))

#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            cmd = "sbatch " +   mdir + '/hf2_aggregate_rasters_submit_' + i + '_' + j + '.sh'
#            print "Submitting Job " + i + " with command: %s" % cmd
#            status, jobnum = commands.getstatusoutput(cmd)
#            if (status == 0):
#                print "Job1 is %s" % jobnum
#            else:
#                print "Error submitting Job " + i

#%%
if runS20 == True:
    ##---- Make hf2_make_convex_hulls_xx.py scripts ----
    repldict = {"scriptname": "hf2_make_convex_hulls", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_make_convex_hulls_submit_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=mch
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/mch_%A_%a.log
    #SBATCH --time=00:03:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=2
    #SBATCH --mem 4000
    #SBATCH --array={njsfupdate}-{njefupdate}
    
    # I think the number of jobs is the number of rows in ppcount.csv minus 1.
    # For big jobs, use 1.5 hours
    
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_make_convex_hulls_{region}_{scenario}.py ${SLURM_ARRAY_TASK_ID}
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Get number of jobarrays
    # and write into the template and then write to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
            odir = '/scratch/pj276/ifl_corridors/' + wdir
            if os.path.exists(odir + '/ppcount.csv'):
                with open(odir + '/ppcount.csv', 'r') as countfile:
                    myjobcount = sum(1 for line in countfile)-1
                # Get number of jobs for update forest cover script using length of file list
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
                with open(mdir + '/hf2_make_convex_hulls_submit_' + i + '_' + j + '.sh','w') as myfile:
                    myfile.write(ttwrite.format(**context))
    
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
#            odir = '/scratch/pj276/ifl_corridors/' + wdir
#            if os.path.exists(odir + '/ppcount.csv'):
#                cmd = "sbatch " +   mdir + '/hf2_make_convex_hulls_submit_' + i + '_' + j + '.sh'
#                print "Submitting Job " + i + " with command: %s" % cmd
#                status, jobnum = commands.getstatusoutput(cmd)
#                if (status == 0):
#                    print "Job1 is %s" % jobnum
#                else:
#                    print "Error submitting Job " + i

#%%
if runS21 == True:
    ##---- Make hf2_get_convex_hull_areas_xx.py scripts ----
    repldict = {"scriptname": "hf2_get_convex_hull_areas", "paramoldname": "param_file", "paramnewname": pbase,
            "syspathold": "/home/pj276/projects/ifl_corridors/code/parameter_files/",
            "syspathnew": "/home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/"}
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            mscript = mainscriptmake(repldict, ssdir, mdir, i=i, j=j)
    
    ##---- Make hf2_get_convex_hull_areas_submit_xx.py scripts ----
    ttwrite = """#!/bin/bash
    #SBATCH --job-name=gcha
    #SBATCH --chdir=/scratch/pj276/ifl_corridors/output
    #SBATCH --output=/scratch/pj276/ifl_corridors/logs/gcha_%A_%a.log
    #SBATCH --time=00:03:00
    #SBATCH --partition=all
    #SBATCH --cpus-per-task=2
    #SBATCH --mem 8000
    #SBATCH --array={njsfupdate}-{njefupdate}
    
    module load anaconda/latest
    
    srun ~/.conda/envs/hifo/bin/python /home/pj276/projects/ifl_corridors/code/scripts/prepnn1_0_2000_3/hf2_get_convex_hull_areas_{region}_{scenario}.py ${SLURM_ARRAY_TASK_ID}
    """
    ttwrite = inspect.cleandoc(ttwrite)
    # Get number of jobarrays
    # and write into the template and then write to file
    for i in ['af','as','sa']:
        for j in ['neg2', 'neg01', 'pos2']:
            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
            odir = '/scratch/pj276/ifl_corridors/' + wdir
            fnlist = []
            if os.path.exists(odir + '/convex_hulls'):
                for dir, subdir, files in os.walk(odir + '/convex_hulls'):
                    for fname in files:
                        if fname.endswith(".shp"):
                            fnlist.extend([fname.split('_')[3]])
                            #fnlist.extend(os.path.join(dir + "/", fname))
                fnlist = list(set(fnlist))
                myjobcount = len(fnlist)
                # Get number of jobs for update forest cover script using length unique grid codes
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
                with open(mdir + '/hf2_get_convex_hull_areas_submit_' + i + '_' + j + '.sh','w') as myfile:
                    myfile.write(ttwrite.format(**context))
    
#    #---- SUBMIT JOBS ----
#    for i in ['af', 'as', 'sa']:
#        for j in ['neg2', 'neg01', 'pos2']:
#            wdir = 'rp5k130k90m_' + j + '_' + i + '_2000_ifl'
#            odir = '/scratch/pj276/ifl_corridors/' + wdir
#            if os.path.exists(odir + '/ppcount.csv'):
#                cmd = "sbatch " +   mdir + '/hf2_get_convex_hull_areas_submit_' + i + '_' + j + '.sh'
#                print "Submitting Job " + i + " with command: %s" % cmd
#                status, jobnum = commands.getstatusoutput(cmd)
#                if (status == 0):
#                    print "Job1 is %s" % jobnum
#                else:
#                    print "Error submitting Job " + i