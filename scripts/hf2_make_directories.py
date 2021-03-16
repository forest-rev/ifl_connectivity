# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 21:50:30 2018
Create directories for corridor mapping
@author: pj276
"""

#%%
import sys, os
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
#import p1_90m_110718_as as p1
#import p1_90m_ifl_2000_as as p1
import param_file as p1

#%%
# Working directory
odir = p1.odir
# Directories to be created
dlist = p1.dlist #['/distfiles/th', '/distfiles/lcp', '/distfiles/ptile', '/cfinal/flux', '/cfinal/ncount', '/cfinal/acd', '/clevel1count', '/clevel1cdist', '/corridors', '/csumfiles1', '/finfiles','/ppc4redo','/convex_hulls','/pointfiles','/areadfs']

# Loop through and create
for i in dlist:
    if not os.path.exists(odir + i):
        os.makedirs(odir + i)
