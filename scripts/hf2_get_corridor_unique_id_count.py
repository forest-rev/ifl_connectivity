# -*- coding: utf-8 -*-
"""
Created on Thu Apr 05 14:05:36 2018
Get corridor count for the summarize corridor function.
@author: pj276
"""
import sys, glob, pandas as pd
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
#import p1_10k_100k_041118 as p1
#import p1_10k_300k_060518 as p1
#import p1_90m_110718_as as p1
#import p1_90m_ifl_2000_sa as p1
import param_file as p1

# Input directory
sdir = p1.sdir
#'/projects/above_gedi/pjantz/ifl_corridors_data/outputs/regular_points/rp10ks1/corridors'

odir = p1.odir
finstatusptn = p1.finstatusptn
ffdir = p1.ffdir

#'corr_neg1_'
# Output text file holding id count 
outidcount = p1.outidcount
# Output file index name
findex = p1.findex
# Output unique corridor pair file
ucpfile = p1.ucpfile

#%%
# List csv files with corridor processing status info
# and concatenate them to a data frame
flist = glob.glob(odir + ffdir + finstatusptn)
flist = [i for i in flist if i.endswith('.csv')]
flist = [pd.read_csv(i) for i in flist]
fdf = pd.concat(flist, ignore_index=True)
# Get only corridors that were successfully processed
fdf = fdf.loc[fdf['status']=='finished']

# Get list of unique patch pairs
fdfgb = fdf.groupby(['gc1', 'gc2'])
fdfgroups = fdfgb.groups.keys()
fdfset = list(set(fdfgroups))

# Write count of patch pairs to file
with open(outidcount, 'w') as afile:
    afile.write(str(len(fdfset)))
    
# Write file list to file so subsequent scripts don't waste time listing them again
fdf.to_csv(findex, index=False)

# Write unique ids to df
adf = pd.DataFrame(columns=['pp1', 'pp2'])
adf['pp1'] = [i[0] for i in fdfset]
adf['pp2'] = [i[1] for i in fdfset]
adf.to_csv(ucpfile, index=False)
