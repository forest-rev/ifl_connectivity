# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 12:40:01 2018
Count the lines in a csv, divide by 10, and write result to file.
This number is used to bound a job array.
@author: pj276
"""

import sys, math
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
#import p1_10k_300k_060518 as p1
#import p1_90m_110718_as as p1
#import p1_90m_ifl_2000_sa as p1
import param_file as p1

countfile = p1.csvout
outcount = p1.outcount

numlines = math.ceil((sum(1 for line in open(countfile))-1)/10.0)

with open(outcount, 'a') as afile:
    afile.write(str(numlines))
