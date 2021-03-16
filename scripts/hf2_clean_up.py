# -*- coding: utf-8 -*-
"""
Created on Thu May  9 15:27:27 2019
Delete areadf_xx.csv and ptnbs_100k_xx.csv files
@author: pj276
"""
# Import modules
import sys, os, shutil

# Path to param files
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
#import p1_90m_110718_as as p1
import param_file as p1

# Set variables
odir = p1.odir
sdir = p1.sdir
areadfs = p1.areadfs

# List area files
a1 = os.listdir(odir + areadfs)
# Get areadf files but omitting areadf.csv
a1 = [i for i in a1 if 'areadf' in i and i != 'areadf.csv']
# Delete files
if len(a1) > 0:
    [os.remove(odir + areadfs + '/' + i) for i in a1]

# List ptnbs files, omitting the important ones
a1 = os.listdir(odir)
a1 = [i for i in a1 if 'ptnbs_100k' in i and i != 'ptnbs_100k.csv' and i != 'ptnbs_100k_all.csv' and i != 'ptnbs_100k_all_orig.csv']
# Delete files
if len(a1) > 0:
    [os.remove(odir + '/' + i) for i in a1]

# Delete 0th level corridors
if os.path.exists(sdir):
    shutil.rmtree(sdir)
