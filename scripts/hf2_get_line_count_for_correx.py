# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 15:51:20 2018
This script gets a line count from a file to use in corridor expander script.
It also sorts rows by distance, smallest to largest.
@author: pj276
"""

import sys, math, os, glob
import pandas as pd
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
#import p1_10k_300k_060518 as p1
#import p1_90m_110718_as as p1
#import p1_90m_ifl_2000_sa as p1
import param_file as p1

countfile = p1.nbcsvall
outcexcount = p1.outcexcount
forfileindex = p1.forfileindex

numlines = math.ceil((sum(1 for line in open(countfile))-1)/250.0)

with open(outcexcount, 'a') as afile:
    afile.write(str(numlines))
    
# Calculate distance between points
# Order rows by patch and point id
f1 = pd.read_csv(countfile)
f1['dist'] = ((f1.cx2-f1.cx1)**2 + (f1.cy2-f1.cy1)**2)**0.5
#f1 = f1.sort_values(by="dist")
f1 = f1.sort_values(by=['gc1','ptid1','gc2','ptid2'])
f1.to_csv(countfile, index=False)

# Create text file for next script if it doesn't exist
if not os.path.exists(forfileindex):
    alist = glob.glob(os.path.dirname(forfileindex) + '/' + '*' + os.path.basename(forfileindex).split('index')[0] + '*')
    with open(forfileindex, 'a') as tfile:
        for k in alist:
            tfile.write(k + '\n')
