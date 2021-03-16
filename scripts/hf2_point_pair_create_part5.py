# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 11:20:48 2018
Concatenate area estimate dataframes
@author: pj276
"""
#%%
import sys, glob
import pandas as pd
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
#import p1_90m_110718_as as p1
import param_file as p1

#%%
# Args from parameter file
odir = p1.odir
adfname_temp = p1.adfname_temp
adfname = p1.adfname

#%%
# List files
flist = glob.glob(adfname_temp + '*')
flist = [i for i in flist if '.csv' in i]
# Empty list
dflist = []
# Read point neighbor files and append to list
for i in (flist):
    adf = pd.read_csv(i)
    dflist.append(adf)

# Concatenate data frames
rdf = pd.concat(dflist, ignore_index=True)

# Get rid of duplicates
rdf = rdf.drop_duplicates()
rdf = rdf.drop(['Unnamed: 0'], axis=1)
rdf = rdf.reset_index()

# Name index
rdf.index.name = 'rid'

# Write to file
rdf.to_csv(adfname)
