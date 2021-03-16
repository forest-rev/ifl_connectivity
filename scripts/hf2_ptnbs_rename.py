# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 16:09:23 2019
Rename neighbors file
@author: pj276
"""

#%%
# Packages
import sys, os
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
import param_file as p1

#%%
# Rename first file
src = p1.srcname1
dst = p1.dstname1

os.rename(src, dst)

# Rename second file
src = p1.srcname2
dst = p1.dstname2

os.rename(src, dst)
