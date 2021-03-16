# -*- coding: utf-8 -*-
"""
Created on Sun Mar 25 22:21:05 2018
This script aggregates csv files from part two of the script to create one
file that gives neighbors for each point.
@author: pj276
"""
#%%
import pandas as pd
import sys, os, glob
import numpy as np
from sklearn import linear_model
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
#import p1_10k_100k_041118 as p1
#import p1_20k_100k_052518 as p1
#import p1_10k_100k_060418 as p1
#import p1_rp_300k_060518 as p1
#import p1_10k_300k_060518 as p1
#import p1_90m_110718_as as p1
#import p1_90m_ifl_2000_sa as p1
import param_file as p1

#%%
#300000 # maximum point distance to consider, in meters
nbcsvall = p1.nbcsvall
#'/projects/above_gedi/pjantz/ifl_corridors_data/outputs/regular_points/rp10ks1/ptnbs_10k_300k_all.csv'                                                                                                             
spat = p1.spat
#'/projects/above_gedi/pjantz/ifl_corridors_data/outputs/regular_points/rp10ks1/ptnbs_10k_300k_*'
# Patch pair count file
ppcount = p1.ppcount + '.csv'
# Patch pair subset file
rdfsub = p1.rdfsub + '.csv'
# Directory
odir = p1.odir

#%%
# Remove temp files from ppc1 script
delist = glob.glob(odir + '/*_temp.csv')
if len(delist) > 0:
    [os.remove(i) for i in delist]

# List files to loop
flist = glob.glob(spat)
flist = [i for i in flist if 'all' not in i]
# Empty list
dflist = []
# Read point neighbor files and append to list
for i in (flist):
    nbcsvout = i
    #nbcsvout = '/projects/above_gedi/pjantz/ifl_corridors_data/outputs/regular_points/rp10ks1/ptnbs_10k_300k_' + str(i) + '.csv'
    adf = pd.read_csv(nbcsvout)
    adf = adf.drop(['rid'],axis=1)
    dflist.append(adf)

# Concatenate data frames
rdf = pd.concat( dflist, ignore_index=True)
# Get rid of duplicates
rdf = rdf.drop_duplicates()
# Remove duplicates across columns
#amask = rdf['uid1'] < rdf['uid2'] # reorder unique ids so they appear in canonical order
#rdf['first'] = rdf['uid1'].where(amask, rdf['uid2'])
#rdf['second'] = rdf['uid2'].where(amask, rdf['uid1'])
#rdf = rdf.drop_duplicates(subset=['first', 'second']) # drop duplicates based on the canonically ordered ids
#rdf = rdf.drop(['uid1', 'uid2', 'first', 'second'], axis=1)

# Name index
rdf.index.name = 'rid'

# Write to file
rdf.to_csv(nbcsvall)

# Get count of the number of connections in each patch pair and write to csv
ppc = rdf.groupby(['gc1','gc2']).size().reset_index().rename(columns={0:'count'})
ppc.index.name = 'rid'
ppc.to_csv(ppcount)

# Subsample large groups down to the threshold
##size = ssize        # sample size
##replace = False  # with replacement
##fn = lambda obj: obj.loc[np.random.choice(obj.index, size, replace),:]
##rdfs = rdf1.groupby(['gc1', 'gc2'], as_index=False).apply(fn)

# Adaptive subsampling
# Set up two vectors to define the relationship between observed
# number of point pairs and the subsample. Convert to log for regression.
# Increase point pairs by doubling
ppairs1 = np.log([100*(2**i) for i in range(0,10)])
# Increase point pairs by 50%
ppairs2 = np.log([100*(1.5**i) for i in range(0,10)])
# Convert to data frames  
df1 = pd.DataFrame(ppairs1, columns=["pcount1"])
df2 = pd.DataFrame(ppairs2, columns=["pcount2"])
# Rename for regression
X = df1
y = df2["pcount2"]
# Log-log linear model
lm = linear_model.LinearRegression()
model = lm.fit(X,y)

# Function to subsample data frame if number of point pairs > 500
# If the new subsampled # of point pairs is less than 500, don't subsample
# This basically means that subsampling occurs when n > 500.
# There's a discontinuity here but I'm not entirely sure how to get rid of it
def ss_x_points(frame,amod):
    x = frame.shape[0]
    newdat = pd.DataFrame.from_dict({'c1':np.log([x])})
    if x > 500 and np.int(np.exp(amod.predict(newdat))) > 500:
        #frame['newcount'] = np.exp(amod.predict(newdat))
        #frame['newcount'] = pd.Series(np.exp(amod.predict(newdat)), index=frame.index)
        np.random.seed(52)
        frame = frame.loc[np.random.choice(frame.index, np.int(np.exp(amod.predict(newdat))), replace = False),:]
        return(frame)
    else:
        #frame['newcount'] = x
        return(frame)

# Create groupby object
dfbg = rdf.groupby(['gc1','gc2'])
# Apply function
newcounts = dfbg.apply(lambda x: ss_x_points(x,lm))
# Create a column from the index
svals = pd.Series(range(0,newcounts.shape[0]))
newcounts['svals'] = svals.values
newcounts = newcounts.set_index('svals')
newcounts.index.name = 'rid'
newcounts.to_csv(rdfsub)

# Create redo directory for the next script
if not os.path.exists(p1.odir + p1.ppc4redodir):
    os.mkdir(p1.odir + p1.ppc4redodir)