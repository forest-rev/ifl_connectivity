# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 14:30:00 2018
This script calculates which points are within a threshold distance. The points
were mapped in part one of the script. It runs in batches of 10.
@author: pj276
"""

#%%
import scipy.spatial as spatial
import pandas as pd
import sys
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
csvout = p1.csvout
#'/projects/above_gedi/pjantz/ifl_corridors_data/outputs/regular_points/rp10ks1/ptnbs_10k.csv'
nbcsvpre = p1.nbcsvpre
#'/projects/above_gedi/pjantz/ifl_corridors_data/outputs/regular_points/rp10ks1/ptnbs_10k_300k_'
tdist = p1.tdist
#300000 # maximum point distance to consider, in meters

# Job array index
j = int(sys.argv[1]) - 1 # subtract for python zero indexing
jstart = 10*j
jend = jstart + 10
nbcsvout = nbcsvpre + str(j) + '.csv'

#%%
# Read source target patches from csv columns to lists
#dat = pandas.read_csv(icsv + str(j) + '.csv')
adf = pd.read_csv(csvout)
gridcodes = adf.gcode.tolist()
xcos = adf.xco.tolist()
ycos = adf.yco.tolist()
ptids = adf.ptid.tolist()
# Create coordinate matrix
cmat = zip(xcos,ycos)
# Create kdtrees
ptree1 = spatial.cKDTree(cmat[jstart:jend])
ptree2 = spatial.cKDTree(cmat)
# Get neighbors within threshold distance (in meters)
# This returns indices relative to ptree2
nbs = ptree1.query_ball_tree(ptree2,r=tdist,p=2)
# Empty lists
nlist = []
#uidlist = []
d = {}

# Change jend if nbs is less than 10. Otherwise it throws an indexing error.
if len(nbs) < 10:
    jend = len(nbs)

# Create dictionary where keys are point-patch ids and elements are x,y coordinates
# Loop through neighbor objects from kdtree and put point pairs
# on a single line. Remove duplicates and intrapatch pairs.
for i,k in enumerate(range(jstart,jend)):
    if len(nbs[i]) > 1:
        for m in nbs[i]:
            if gridcodes[k] != gridcodes[m]:
                d[str(gridcodes[k]) + '_' + str(ptids[k])] = (xcos[k],ycos[k])
                d[str(gridcodes[m]) + '_' + str(ptids[m])] = (xcos[m],ycos[m])
                nlist.append([gridcodes[k], ptids[k], xcos[k], ycos[k], gridcodes[m], ptids[m], xcos[m], ycos[m]])
#            uid1 = [str(gridcodes[k]) + '_' + str(ptids[k]) + '_' + str(gridcodes[m]) + '_' + str(ptids[m])]
#            uid2 = [str(gridcodes[m]) + '_' + str(ptids[m]) + '_' + str(gridcodes[k]) + '_' + str(ptids[k])]
#            if gridcodes[k] != gridcodes[m] and uid1 not in uidlist and uid2 not in uidlist:
#                nlist.append([gridcodes[k], ptids[k], xcos[k], ycos[k], gridcodes[m], ptids[m], xcos[m], ycos[m]])
#                uidlist.append(uid1)
#                uidlist.append(uid2)

# Dump to data frame
nbdf = pd.DataFrame(nlist)
try:
    # Name columns
    nbdf.columns = ['gc1', 'ptid1', 'cx1', 'cy1', 'gc2', 'ptid2', 'cx2', 'cy2']
except ValueError:
    with open(p1.odir + '/ppc2errors_' + str(j) + '.txt', 'a') as tf:
        tf.write('Value error in file ' + str(j) + '\n')
    sys.exit(1)
# Combine patch and point ids for each patch-point pair
nbdf['id1'] = nbdf['gc1'].map(str) + '_' + nbdf['ptid1'].map(str)
nbdf['id2'] = nbdf['gc2'].map(str) + '_' + nbdf['ptid2'].map(str)
# Calculate mins and maxes across columns to derive a consistent ordering
# of patch-point pairs. E.g. given pp1: 1876_30 and 2030_5; pp2: 2030_5 and 1876_30
# pp2 would be reordered to 1876_30 and 2030_5
nbdf['id_min'] = nbdf.loc[:, ['id1', 'id2']].min(axis=1)
nbdf['id_max'] = nbdf.loc[:, ['id1', 'id2']].max(axis=1)
# Split patch-point pairs back to individual columns
nbdf['gc1'], nbdf['ptid1'] = nbdf['id_min'].str.split('_', 1).str
nbdf['gc2'], nbdf['ptid2'] = nbdf['id_max'].str.split('_', 1).str
# Drop coordinate columns as their relationship with patch-point IDs has been broken
nbdf = nbdf.drop(['cx1', 'cx2', 'cy1', 'cy2', 'id1', 'id2'], 1)

# Loop through list of first point-patch ids and append associated coords to list
xlist1 = []
ylist1 = []
for k in list(nbdf['id_min']):
    if k in d:
        xlist1.append(d[k][0])
        ylist1.append(d[k][1])

# Loop through list of second point-patch ids and append associated coords to list
xlist2 = []
ylist2 = []
for k in list(nbdf['id_max']):
    if k in d:
        xlist2.append(d[k][0])
        ylist2.append(d[k][1])

# Add coordinates back to data frame
nbdf['cx1'] = xlist1
nbdf['cy1'] = ylist1
nbdf['cx2'] = xlist2
nbdf['cy2'] = ylist2

# Drop extra id columns
nbdf = nbdf.drop(['id_min', 'id_max'], 1)
    
# Remove duplicates
nbdf = nbdf.drop_duplicates()
# Name index
nbdf.index.name = 'rid'
# Write to file
nbdf.to_csv(nbcsvout)


