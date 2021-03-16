# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 13:37:44 2019
Summarize forest cover in IFL corridors
@author: pj276
"""

#%%
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib as plt
import glob
from decimal import *

%matplotlib qt

#%%
# fc by region
# List files
flist = glob.glob('G:/My Drive/Projects/ifl_corridors/statsfiles/fc_region_stats/*')
flist = [pd.read_csv(f) for f in flist]
# Concatenate files
dfs = pd.concat(flist)
# Split by year and sort
df1 = dfs[dfs['year'] == 2000]
df2 = dfs[dfs['year'] == 2013]
df1 = df1.sort_values(by=['region', 'scenario'])
df2 = df2.sort_values(by=['region', 'scenario'])
# Calculate absolute difference over time
absdiff = df1.meanc - df2.meanc
# Percent difference
pctdiff = (df2.meanc - df1.meanc)/df1.meanc*100

# Get summary of regional forest cover
regfc = pd.concat([df1.meanc,df1.sdc,df2.meanc,df2.sdc,absdiff, pctdiff, df1.region, df1.scenario],
                  keys=['fc2000','fcsd2000','fc2013','fcsd2013','absdiff','pdiff','region','scenario'],
                  axis=1, ignore_index=False)


#%%
# fc by country
# List files
flist = glob.glob('G:/My Drive/Projects/ifl_corridors/statsfiles/gadm_fc/*')
flist = [pd.read_csv(f) for f in flist]
# Concatenate files
dfs = pd.concat(flist)
# Remove rows with no values
dfs = dfs.loc[dfs['fcave'] != '--']
# Convert columns to float
dfs = dfs.astype({'fcave': 'float32', 'fcstd': 'float32'})
# Rename ucode column to country
dfs.rename(columns={'ucode': 'country'}, inplace=True)
# Subset out neg01 scenario
dfs = dfs.loc[dfs['scenario']=='neg01']
# Change 2013 to 2012 in year column
dfs.loc[dfs['year']==2013,'year'] = 2012

# Split by year and sort
df1 = dfs[dfs['year'] == 2000]
df2 = dfs[dfs['year'] == 2012]

df1 = df1.sort_values(by=['region', 'scenario', 'country'])
df2 = df2.sort_values(by=['region', 'scenario', 'country'])
df1.reset_index(inplace = True)
df2.reset_index(inplace = True)
# Calculate difference over time
absdiff = df1['fcave'] - df2['fcave']
pctdiff = (df1.fcave - df2.fcave)/df1.fcave*100

# Get list of countries sorted by change in forest cover in corridors
csort = [x for _,x in sorted(zip(absdiff,df1.country.tolist()), reverse=True)]

#ax = sns.catplot(x='year',y='fcave',hue='country',kind='point',col='scenario',data=dfs)
#ax = sns.catplot(x='country',y='fcave',hue='year',kind='bar',col='scenario',data=dfs)
#ax = sns.catplot(x='country',y='fcave',hue='year',kind='bar',data=dfs, order=csort, palette=("Paired"))
sns.set_context('poster',rc={"font.size":40,"axes.titlesize":40,"axes.labelsize":40})
ax = sns.catplot(x='country',y='fcave',hue='year',kind='bar',data=dfs, order=csort, palette=(['#6b8ba4','#fac205']))
ax.set_xticklabels(rotation=45, horizontalalignment='right')
ax.set(xlabel='Country', ylabel='% Tree Cover')




