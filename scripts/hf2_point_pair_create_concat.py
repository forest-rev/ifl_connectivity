# -*- coding: utf-8 -*-
"""
Created on Sun Aug 12 13:08:14 2018
This script concatentates csv files of point locations.
@author: pj276
"""
import pandas as pd
import sys, glob
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
#import p1_90m_110718_as as p1
#import p1_90m_ifl_2000_sa as p1
import param_file as p1

# Set up inputs
odir = p1.odir
csvout = p1.csvout
fpre = p1.fpre # ptnbs_100k

# List csvs, join, and save.
flist = glob.glob(odir + '/' + fpre + '_*')
flist = [i for i in flist if '.csv' in i]
elist = []
for z in flist:
    elist.append(pd.read_csv(z))
elist = pd.concat(elist)
elist = elist.drop(['rid'], axis=1)
elist.reset_index(drop=True, inplace=True)
elist.index.name = 'rid'
elist.to_csv(csvout)
