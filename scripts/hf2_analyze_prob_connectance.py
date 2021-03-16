# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 13:08:07 2019
Get flux change stats, connectance change stats, area/matrix contribution to flux
@author: pj276
"""

#%%
import pandas as pd
import numpy as np
import glob

#%matplotlib qt

#%%
path = 'G:/My Drive/Projects/ifl_corridors/rp5k130k25k/dstats_all'

#-----------------------------------------------------------------------------------------
## 1 ##
# Calculate overall change in connectivity and by region
#%%
## pos2 ##
subdirs = ['dstatsl2_rp5k130k90m_pos2_af_2000_ifl','dstatsl2_rp5k130k90m_pos2_as_2000_ifl','dstatsl2_rp5k130k90m_pos2_sa_2000_ifl',
           'dstatsl2_rp5k130k90m_pos2_af_2013_ifl','dstatsl2_rp5k130k90m_pos2_as_2013_ifl','dstatsl2_rp5k130k90m_pos2_sa_2013_ifl']
dfsall = []
for i in subdirs:
    flist = glob.glob(path + '/' + i + '/*' )
    dfs = [pd.read_csv(j) for j in flist]
    dfs = pd.concat(dfs)
    if '2000' in i:
        dfs['year'] = '2000'
    if '2013' in i:
        dfs['year'] = '2013'
    if 'af' in i:
        dfs['region'] = 'af'
    if 'as' in i:
        dfs['region'] = 'as'
    if 'sa' in i:
        dfs['region'] = 'sa'
    dfsall.append(dfs)
dfsall = pd.concat(dfsall)

# Remove corridor pairs with no flux
dfsall = dfsall.loc[dfsall['p10probave']!=0]

# Get average probability of flux change
pos2fsum = dfsall.groupby(['region','year'])['p10fluxave'].sum()
# Africa
pos2affch = (pos2fsum['af','2013'] - pos2fsum['af','2000'])/pos2fsum['af','2000']*100
# Asia
pos2asfch = (pos2fsum['as','2013'] - pos2fsum['as','2000'])/pos2fsum['as','2000']*100
# South America
pos2safch = (pos2fsum['sa','2013'] - pos2fsum['sa','2000'])/pos2fsum['sa','2000']*100

# Pan-tropical
pos2fsum = dfsall.groupby(['year'])['p10fluxave'].sum()
pos2ptsum = (pos2fsum['2013'] - pos2fsum['2000'])/pos2fsum['2000']*100

## neg01 ##
subdirs = ['dstatsl2_rp5k130k90m_neg01_af_2000_ifl','dstatsl2_rp5k130k90m_neg01_as_2000_ifl','dstatsl2_rp5k130k90m_neg01_sa_2000_ifl',
           'dstatsl2_rp5k130k90m_neg01_af_2013_ifl','dstatsl2_rp5k130k90m_neg01_as_2013_ifl','dstatsl2_rp5k130k90m_neg01_sa_2013_ifl']
dfsall = []
for i in subdirs:
    flist = glob.glob(path + '/' + i + '/*' )
    dfs = [pd.read_csv(j) for j in flist]
    dfs = pd.concat(dfs)
    if '2000' in i:
        dfs['year'] = '2000'
    if '2013' in i:
        dfs['year'] = '2013'
    if 'af' in i:
        dfs['region'] = 'af'
    if 'as' in i:
        dfs['region'] = 'as'
    if 'sa' in i:
        dfs['region'] = 'sa'
    dfsall.append(dfs)
dfsall = pd.concat(dfsall)

# Remove corridor pairs with no flux
dfsall = dfsall.loc[dfsall['p10probave']!=0]

# Get average probability of flux change
neg01fsum = dfsall.groupby(['region','year'])['p10fluxave'].sum()
# Africa
neg01affch = (neg01fsum['af','2013'] - neg01fsum['af','2000'])/neg01fsum['af','2000']*100
# Asia
neg01asfch = (neg01fsum['as','2013'] - neg01fsum['as','2000'])/neg01fsum['as','2000']*100
# South America
neg01safch = (neg01fsum['sa','2013'] - neg01fsum['sa','2000'])/neg01fsum['sa','2000']*100

# Pan-tropical
neg01fsum = dfsall.groupby(['year'])['p10fluxave'].sum()
neg01ptsum = (neg01fsum['2013'] - neg01fsum['2000'])/neg01fsum['2000']*100


## neg2 ##
subdirs = ['dstatsl2_rp5k130k90m_neg2_af_2000_ifl','dstatsl2_rp5k130k90m_neg2_as_2000_ifl','dstatsl2_rp5k130k90m_neg2_sa_2000_ifl',
           'dstatsl2_rp5k130k90m_neg2_af_2013_ifl','dstatsl2_rp5k130k90m_neg2_as_2013_ifl','dstatsl2_rp5k130k90m_neg2_sa_2013_ifl']
dfsall = []
for i in subdirs:
    flist = glob.glob(path + '/' + i + '/*' )
    dfs = [pd.read_csv(j) for j in flist]
    dfs = pd.concat(dfs)
    if '2000' in i:
        dfs['year'] = '2000'
    if '2013' in i:
        dfs['year'] = '2013'
    if 'af' in i:
        dfs['region'] = 'af'
    if 'as' in i:
        dfs['region'] = 'as'
    if 'sa' in i:
        dfs['region'] = 'sa'
    dfsall.append(dfs)
dfsall = pd.concat(dfsall)

# Remove corridor pairs with no flux
dfsall = dfsall.loc[dfsall['p10probave']!=0]

# Get average probability of connectance change
neg2fsum = dfsall.groupby(['region','year'])['p10fluxave'].sum()
# Africa
neg2affch = (neg2fsum['af','2013'] - neg2fsum['af','2000'])/neg2fsum['af','2000']*100
# Asia
neg2asfch = (neg2fsum['as','2013'] - neg2fsum['as','2000'])/neg2fsum['as','2000']*100
# South America
neg2safch = (neg2fsum['sa','2013'] - neg2fsum['sa','2000'])/neg2fsum['sa','2000']*100

# Pan-tropical
neg2fsum = dfsall.groupby(['year'])['p10fluxave'].sum()
neg2ptsum = (neg2fsum['2013'] - neg2fsum['2000'])/neg2fsum['2000']*100

## Summaries ##
print 'pos2 percent change in connectivity'
print pos2ptsum

print 'neg01 percent change in connectivity'
print neg01ptsum

print 'neg2 percent change in connectivity'
print neg2ptsum

#-----------------------------------------------------------------------------------------
## 2 ##
# Summarize probability of connectance and euclidean distance between patches in each region
#%%

## pos2 ##
subdirs = ['dstatsl1_rp5k130k90m_pos2_af_2000_ifl','dstatsl1_rp5k130k90m_pos2_as_2000_ifl','dstatsl1_rp5k130k90m_pos2_sa_2000_ifl',
           'dstatsl1_rp5k130k90m_pos2_af_2013_ifl','dstatsl1_rp5k130k90m_pos2_as_2013_ifl','dstatsl1_rp5k130k90m_pos2_sa_2013_ifl']
dfsall = []
for i in subdirs:
    flist = glob.glob(path + '/' + i + '/*' )
    dfs = [pd.read_csv(j) for j in flist]
    dfs = pd.concat(dfs)
    if '2000' in i:
        dfs['year'] = '2000'
    if '2013' in i:
        dfs['year'] = '2013'
    if 'af' in i:
        dfs['region'] = 'af'
    if 'as' in i:
        dfs['region'] = 'as'
    if 'sa' in i:
        dfs['region'] = 'sa'
    dfsall.append(dfs)
dfsall = pd.concat(dfsall)

# Get average probability of connectance
pos2aprob = dfsall.groupby(['region','year'])['apd'].mean()
print 'average probability of connectance for pos2'
print pos2aprob
# Get average euclidean distance
pos2adist = dfsall.groupby(['region','year'])['dist'].mean()
print 'average euclidean distance for pos2'
print pos2adist


# Clean up
del dfsall

## neg01 ##
subdirs = ['dstatsl1_rp5k130k90m_neg01_af_2000_ifl','dstatsl1_rp5k130k90m_neg01_as_2000_ifl','dstatsl1_rp5k130k90m_neg01_sa_2000_ifl',
           'dstatsl1_rp5k130k90m_neg01_af_2013_ifl','dstatsl1_rp5k130k90m_neg01_as_2013_ifl','dstatsl1_rp5k130k90m_neg01_sa_2013_ifl']
dfsall = []
for i in subdirs:
    flist = glob.glob(path + '/' + i + '/*' )
    dfs = [pd.read_csv(j) for j in flist]
    dfs = pd.concat(dfs)
    if '2000' in i:
        dfs['year'] = '2000'
    if '2013' in i:
        dfs['year'] = '2013'
    if 'af' in i:
        dfs['region'] = 'af'
    if 'as' in i:
        dfs['region'] = 'as'
    if 'sa' in i:
        dfs['region'] = 'sa'
    dfsall.append(dfs)
dfsall = pd.concat(dfsall)

# Get average probability of connectance
neg01aprob = dfsall.groupby(['region','year'])['apd'].mean()
print 'average probability of connectance for neg01'
print neg01aprob
# Get average euclidean distance
neg01adist = dfsall.groupby(['region','year'])['dist'].mean()
print 'average euclidean distance for neg01'
print neg01adist

# Clean up
del dfsall

## neg2 ##
subdirs = ['dstatsl1_rp5k130k90m_neg2_af_2000_ifl','dstatsl1_rp5k130k90m_neg2_as_2000_ifl','dstatsl1_rp5k130k90m_neg2_sa_2000_ifl',
           'dstatsl1_rp5k130k90m_neg2_af_2013_ifl','dstatsl1_rp5k130k90m_neg2_as_2013_ifl','dstatsl1_rp5k130k90m_neg2_sa_2013_ifl']
dfsall = []
for i in subdirs:
    flist = glob.glob(path + '/' + i + '/*' )
    dfs = [pd.read_csv(j) for j in flist]
    dfs = pd.concat(dfs)
    if '2000' in i:
        dfs['year'] = '2000'
    if '2013' in i:
        dfs['year'] = '2013'
    if 'af' in i:
        dfs['region'] = 'af'
    if 'as' in i:
        dfs['region'] = 'as'
    if 'sa' in i:
        dfs['region'] = 'sa'
    dfsall.append(dfs)
dfsall = pd.concat(dfsall)

# Get average probability of connectance
neg2aprob = dfsall.groupby(['region','year'])['apd'].mean()
print 'average probability of connectance for neg2'
print neg2aprob
# Get average euclidean distance
neg2adist = dfsall.groupby(['region','year'])['dist'].mean()
print 'average euclidean distance for neg2'
print neg2adist

# Clean up
del dfsall

#-----------------------------------------------------------------------------------------
## 3 ##
# Calculate share of flux change due to change in area or change in probability of connectance
#%%
## pos2 ##
subdirs = ['dstatsl2_rp5k130k90m_pos2_af_2000_ifl','dstatsl2_rp5k130k90m_pos2_as_2000_ifl','dstatsl2_rp5k130k90m_pos2_sa_2000_ifl',
           'dstatsl2_rp5k130k90m_pos2_af_2013_ifl','dstatsl2_rp5k130k90m_pos2_as_2013_ifl','dstatsl2_rp5k130k90m_pos2_sa_2013_ifl']
dfsall = []
for i in subdirs:
    flist = glob.glob(path + '/' + i + '/*' )
    dfs = [pd.read_csv(j) for j in flist]
    dfs = pd.concat(dfs)
    if '2000' in i:
        dfs['year'] = '2000'
    if '2013' in i:
        dfs['year'] = '2013'
    if 'af' in i:
        dfs['region'] = 'af'
    if 'as' in i:
        dfs['region'] = 'as'
    if 'sa' in i:
        dfs['region'] = 'sa'
    dfsall.append(dfs)
dfsall = pd.concat(dfsall)

# Get corridor heterogeneity var 1
pos2chet1 = dfsall.groupby(['region','year'])['cprobVar1'].mean()
print 'average corridor heterogeneity var1 for pos2'
print pos2chet1
# Get corridor heterogeneity var 2
pos2chet2 = dfsall.groupby(['region','year'])['cprobVar2'].mean()
print 'average corridor heterogeneity var2 for pos2'
print pos2chet2

# Get a data frame for each year
df1 = dfsall[dfsall['year'] == '2000']
df2 = dfsall[dfsall['year'] == '2013']

# Merge them
gnew = pd.merge(df1, df2, how='inner', on=None, left_on=['gc1','gc2','region'], right_on=['gc1','gc2','region'],
         left_index=False, right_index=False, sort=False,
         suffixes=('_1', '_2'), copy=True, indicator=False,
         validate=None)

# Get observed flux difference
gnew['obsdiff'] = gnew['p10fluxave_2']-gnew['p10fluxave_1']
# Get flux difference if probability of connectance did not change
gnew['nchdiff'] = (gnew['gc1Area_2']*gnew['gc2Area_2']*gnew['p10probave_1'])-gnew['p10fluxave_1']

# Get contribution due to area
gnew['cont_area'] = gnew['nchdiff']/gnew['obsdiff']*100
# Get contribution due to connectance
gnew['cont_pc'] = 100-gnew['cont_area']

oa1 = gnew.groupby(['region'])['obsdiff'].sum()
oa2 = gnew.groupby(['region'])['nchdiff'].sum()

print 'contribution due to area change pos2'
print oa2/oa1*100
print 'contribution due to connectivity change pos2'
print 100-(oa2/oa1*100)

print 'Number of pairs where connectivity increased pos2'
print np.sum(gnew['obsdiff'] > 0)

print 'Number of pairs where connectivity decreased pos2'
print np.sum(gnew['obsdiff'] < 0)

print 'Number of pairs with no change pos2'
print np.sum(gnew['obsdiff'] == 0)

print 'Fraction of pairs where connectivity increased pos2'
print np.float(np.sum(gnew['obsdiff'] > 0)) / (np.sum(gnew['obsdiff'] > 0) + np.sum(gnew['obsdiff'] < 0) + np.sum(gnew['obsdiff'] == 0))


## neg01 ##
subdirs = ['dstatsl2_rp5k130k90m_neg01_af_2000_ifl','dstatsl2_rp5k130k90m_neg01_as_2000_ifl','dstatsl2_rp5k130k90m_neg01_sa_2000_ifl',
           'dstatsl2_rp5k130k90m_neg01_af_2013_ifl','dstatsl2_rp5k130k90m_neg01_as_2013_ifl','dstatsl2_rp5k130k90m_neg01_sa_2013_ifl']
dfsall = []
for i in subdirs:
    flist = glob.glob(path + '/' + i + '/*' )
    dfs = [pd.read_csv(j) for j in flist]
    dfs = pd.concat(dfs)
    if '2000' in i:
        dfs['year'] = '2000'
    if '2013' in i:
        dfs['year'] = '2013'
    if 'af' in i:
        dfs['region'] = 'af'
    if 'as' in i:
        dfs['region'] = 'as'
    if 'sa' in i:
        dfs['region'] = 'sa'
    dfsall.append(dfs)
dfsall = pd.concat(dfsall)

# Get corridor heterogeneity var 1
neg01chet1 = dfsall.groupby(['region','year'])['cprobVar1'].mean()
print 'average corridor heterogeneity var1 for neg01'
print neg01chet1
# Get corridor heterogeneity var 2
neg01chet2 = dfsall.groupby(['region','year'])['cprobVar2'].mean()
print 'average corridor heterogeneity var2 for neg01'
print neg01chet2

# Get a data frame for each year
df1 = dfsall[dfsall['year'] == '2000']
df2 = dfsall[dfsall['year'] == '2013']

# Merge them
gnew = pd.merge(df1, df2, how='inner', on=None, left_on=['gc1','gc2','region'], right_on=['gc1','gc2','region'],
         left_index=False, right_index=False, sort=False,
         suffixes=('_1', '_2'), copy=True, indicator=False,
         validate=None)

# Get observed flux difference
gnew['obsdiff'] = gnew['p10fluxave_2']-gnew['p10fluxave_1']
# Get flux difference if probability of connectance did not change
gnew['nchdiff'] = (gnew['gc1Area_2']*gnew['gc2Area_2']*gnew['p10probave_1'])-gnew['p10fluxave_1']

# Get contribution due to area
gnew['cont_area'] = gnew['nchdiff']/gnew['obsdiff']*100
# Get contribution due to connectance
gnew['cont_pc'] = 100-gnew['cont_area']

oa1 = gnew.groupby(['region'])['obsdiff'].sum()
oa2 = gnew.groupby(['region'])['nchdiff'].sum()

print 'contribution due to area change neg01'
print oa2/oa1*100
print 'contribution due to connectivity change neg01'
print 100-(oa2/oa1*100)

print 'Number of pairs where connectivity increased neg01'
print np.sum(gnew['obsdiff'] > 0)

print 'Number of pairs where connectivity decreased neg01'
print np.sum(gnew['obsdiff'] < 0)

print 'Number of pairs with no change neg01'
print np.sum(gnew['obsdiff'] == 0)

print 'Fraction of pairs where connectivity increased neg01'
print np.float(np.sum(gnew['obsdiff'] > 0)) / (np.sum(gnew['obsdiff'] > 0) + np.sum(gnew['obsdiff'] < 0) + np.sum(gnew['obsdiff'] == 0))


## neg2 ##
subdirs = ['dstatsl2_rp5k130k90m_neg2_af_2000_ifl','dstatsl2_rp5k130k90m_neg2_as_2000_ifl','dstatsl2_rp5k130k90m_neg2_sa_2000_ifl',
           'dstatsl2_rp5k130k90m_neg2_af_2013_ifl','dstatsl2_rp5k130k90m_neg2_as_2013_ifl','dstatsl2_rp5k130k90m_neg2_sa_2013_ifl']
dfsall = []
for i in subdirs:
    flist = glob.glob(path + '/' + i + '/*' )
    dfs = [pd.read_csv(j) for j in flist]
    dfs = pd.concat(dfs)
    if '2000' in i:
        dfs['year'] = '2000'
    if '2013' in i:
        dfs['year'] = '2013'
    if 'af' in i:
        dfs['region'] = 'af'
    if 'as' in i:
        dfs['region'] = 'as'
    if 'sa' in i:
        dfs['region'] = 'sa'
    dfsall.append(dfs)
dfsall = pd.concat(dfsall)

# Get corridor heterogeneity var 1
neg2chet1 = dfsall.groupby(['region','year'])['cprobVar1'].mean()
print 'average corridor heterogeneity var1 for neg2'
print neg2chet1
# Get corridor heterogeneity var 2
neg2chet2 = dfsall.groupby(['region','year'])['cprobVar2'].mean()
print 'average corridor heterogeneity var2 for neg2'
print neg2chet2

# Get a data frame for each year
df1 = dfsall[dfsall['year'] == '2000']
df2 = dfsall[dfsall['year'] == '2013']

# Merge them
gnew = pd.merge(df1, df2, how='inner', on=None, left_on=['gc1','gc2','region'], right_on=['gc1','gc2','region'],
         left_index=False, right_index=False, sort=False,
         suffixes=('_1', '_2'), copy=True, indicator=False,
         validate=None)

# Get observed flux difference
gnew['obsdiff'] = gnew['p10fluxave_2']-gnew['p10fluxave_1']
# Get flux difference if probability of connectance did not change
gnew['nchdiff'] = (gnew['gc1Area_2']*gnew['gc2Area_2']*gnew['p10probave_1'])-gnew['p10fluxave_1']

# Get contribution due to area
gnew['cont_area'] = gnew['nchdiff']/gnew['obsdiff']*100
# Get contribution due to connectance
gnew['cont_pc'] = 100-gnew['cont_area']

oa1 = gnew.groupby(['region'])['obsdiff'].sum()
oa2 = gnew.groupby(['region'])['nchdiff'].sum()

print 'contribution due to area change neg2'
print oa2/oa1*100
print 'contribution due to connectivity change neg2'
print 100-(oa2/oa1*100)

print 'Number of pairs where connectivity increased neg2'
print np.sum(gnew['obsdiff'] > 0)

print 'Number of pairs where connectivity decreased neg2'
print np.sum(gnew['obsdiff'] < 0)

print 'Number of pairs with no change neg2'
print np.sum(gnew['obsdiff'] == 0)

print 'Fraction of pairs where connectivity increased neg2'
print np.float(np.sum(gnew['obsdiff'] > 0)) / (np.sum(gnew['obsdiff'] > 0) + np.sum(gnew['obsdiff'] < 0) + np.sum(gnew['obsdiff'] == 0))

0.045658/0.002543