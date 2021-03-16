# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 21:38:52 2020

@author: pj276
"""

#%%
import pandas as pd
import seaborn as sns
import numpy as np
import os, glob
import matplotlib
import matplotlib.pyplot as plt

#%%
# Input dir for csvs
idir = 'G:/My Drive/Projects/ifl_corridors/statsfiles/lybytc2000_fluxes/'

%matplotlib auto

#------------------------------------------------------------
# BY COUNTRY
#------------------------------------------------------------
#%%
flist = glob.glob(os.path.join(idir, '*lybytc2000_flux_sums.csv'))
flist = [pd.read_csv(f) for f in flist]
# Concatenate files
dfs = pd.concat(flist)
# Delete extraneous columns
dfs = dfs.drop(columns=['Unnamed: 0', 'Unnamed: 0.1'])
# Rename ucode to country
dfs = dfs.rename(columns={"ucode": "country"})
# Remove rows with no flux in 2000
dfs = dfs[dfs.flux2000 > 0]

flist = glob.glob(os.path.join(idir, '*lybytc2000noifl_flux_sums.csv'))
flist = [pd.read_csv(f) for f in flist]
# Concatenate files
dfs2 = pd.concat(flist)
# Delete extraneous columns
dfs2 = dfs2.drop(columns=['Unnamed: 0', 'Unnamed: 0.1'])
# Rename ucode to country
dfs2 = dfs2.rename(columns={"ucode": "country"})

dfs = pd.merge(dfs, dfs2, how='inner', on=None, left_on=['region','scenario','country','year'], right_on=['region','scenario','country','year'],
         left_index=False, right_index=False, sort=False,
         suffixes=('_1', '_2'), copy=True, indicator=False,
         validate=None)
#!! temporary transformation
dfs.tcacnoifl = dfs.tcacnoifl*100/8100

#%%
# Take floss column that is a string representation of a list
# and explode to multiple columns
floss = dfs['floss'].apply(lambda x: pd.Series(map(float, pd.Series(x).to_list()[0].strip('][').split(', '))))
# Rename columns
floss = floss.rename(columns = lambda x : 'floss_' + str(x))
# Add back to main data frame
dfs = pd.concat([dfs[:], floss[:]], axis=1)
# Remove columns
dfs = dfs.drop(columns=['floss'])

#%%
# Take tcloss column that is a string representation of a list
# and explode to multiple columns
tcloss = dfs['tcloss'].apply(lambda x: pd.Series(map(float, pd.Series(x).to_list()[0].strip('][').split(', '))))
# Rename columns
tcloss = tcloss.rename(columns = lambda x : 'tcloss_' + str(x))
# Add back to main data frame
dfs = pd.concat([dfs[:], tcloss[:]], axis=1)
# Remove columns
dfs = dfs.drop(columns=['tcloss'])

#%%
# Take tclossc column that is a string representation of a list
# and explode to multiple columns
tclossc = dfs['tclossc'].apply(lambda x: pd.Series(map(float, pd.Series(x).to_list()[0].strip('][').split(', '))))
# Rename columns
tclossc = tclossc.rename(columns = lambda x : 'tclossc_' + str(x))
# Add back to main data frame
dfs = pd.concat([dfs[:], tclossc[:]], axis=1)
# Remove columns
dfs = dfs.drop(columns=['tclossc'])

#mycols = ['floss_0', 'floss_1', 'floss_2', 'floss_3', 'floss_4', 'floss_5',
#          'floss_6', 'floss_7', 'floss_8', 'floss_9', 'floss_10', 'floss_11',
#          'floss_12', 'floss_13', 'floss_14', 'floss_15', 'floss_16',
#          'floss_17', 'floss_18']
#
#dfs['fsum'] = dfs.loc[:, mycols].sum(axis=1)
#dfss = dfs.loc[:, mycols].div(dfs['flux2000'], axis=0)
#dfss = dfss*100
#dfss['country'] = dfs.country
#
#mycols = ['tcloss_0', 'tcloss_1', 'tcloss_2', 'tcloss_3', 'tcloss_4', 'tcloss_5',
#          'tcloss_6', 'tcloss_7', 'tcloss_8', 'tcloss_9', 'tcloss_10', 'tcloss_11',
#          'tcloss_12', 'tcloss_13', 'tcloss_14', 'tcloss_15', 'tcloss_16',
#          'tcloss_17', 'tcloss_18']
#
#dfs['tcsum'] = dfs.loc[:, mycols].sum(axis=1)
#dfss = dfs.loc[:, mycols].div(dfs['tca2000'], axis=0)
#dfss = dfss*100
#dfss['country'] = dfs.country

#%%
# Make new columns to hold the ratio of flux lost to corridor tree cover loss
cols = ['ftcratio_{:01d}'.format(i) for i in np.arange(0, 18+1)]
# Calcualte ratio
dfs[cols] = dfs.filter(regex='^floss_').div(dfs.filter(regex='^tclossc_').values)

#%%
# Melt so annual metrics plot as a line for each country
mdat = dfs.melt(['country','scenario', 'region', 'year'], var_name='cols',  value_name='vals')
mdat = mdat.astype({'vals': float})

# Get indices with flux loss and tree cover loss in corridor and the ratio between them 
idx = list(set([i for i in mdat.cols if any(['floss' in i, 'tclossc' in i, 'ftcratio' in i])]))
# Remove floss_0 and tclossc_0
idx = [i for i in idx if all([i != 'floss_0', i != 'tclossc_0', i != 'ftcratio_0'])]
# Get rows of flux and tree cover loss in corridors and the ratio between them
df1 = mdat[mdat['cols'].isin(idx)]
# Get neg01 scenario
df1 = df1[df1['scenario'].isin(['neg01'])]
# Create new numeric column for flux, tree cover loss, and the ratio between them
df1 = df1.assign(yrs = np.tile(np.repeat([range(1,19,1)], np.sum(df1.cols=='floss_1')),3))

# Create label column for flux and tree cover loss in corridors and the ratio between them
def f(row):
    if 'floss' in row['cols']:
        val = 'floss'
    elif 'tclossc' in row['cols']:
        val = 'tclossc'
    elif 'ftcratio' in row['cols']:
        val = 'ftcratio'
    else:
        val = 'notype'
    return val
df1['losstype'] = df1.apply(f, axis=1)

#palette = dict(zip(dots.coherence.unique(),
#                   sns.color_palette("rocket_r", 6)))
#g = sns.relplot(x='yrs', y='vals', hue='country', data=df1, kind='line', col='region')

g = sns.relplot(x='yrs', y='vals', hue='country', data=df1, kind='line', col='losstype', row='region', facet_kws={'sharey': False, 'sharex': False})

g.axes[2,0].set_ylim(0,2000000)
g.axes[0,1].set_xlim(-40,10)


#%%
#------------------------------------------------------------
# BY REGION
#------------------------------------------------------------

#%%
flist = glob.glob(os.path.join(idir, '*lybytc2000_flux_sums.csv'))
flist = [pd.read_csv(f) for f in flist]
# Concatenate files
dfs = pd.concat(flist)
# Delete extraneous columns
dfs = dfs.drop(columns=['Unnamed: 0', 'Unnamed: 0.1'])
# Rename ucode to country
dfs = dfs.rename(columns={"ucode": "country"})
# Remove rows with no flux in 2000
dfs = dfs[dfs.flux2000 > 0]

flist = glob.glob(os.path.join(idir, '*lybytc2000noifl_flux_sums.csv'))
flist = [pd.read_csv(f) for f in flist]
# Concatenate files
dfs2 = pd.concat(flist)
# Delete extraneous columns
dfs2 = dfs2.drop(columns=['Unnamed: 0', 'Unnamed: 0.1'])
# Rename ucode to country
dfs2 = dfs2.rename(columns={"ucode": "country"})

dfs = pd.merge(dfs, dfs2, how='inner', on=None, left_on=['region','scenario','country','year'], right_on=['region','scenario','country','year'],
         left_index=False, right_index=False, sort=False,
         suffixes=('_1', '_2'), copy=True, indicator=False,
         validate=None)
#!! temporary transformation
dfs.tcacnoifl = dfs.tcacnoifl*100/8100

#df2 = df1.groupby(['region', 'yrs', 'losstype'], as_index=False).sum()

#%%
# Take floss column that is a string representation of a list
# and explode to multiple columns
floss = dfs['floss'].apply(lambda x: pd.Series(map(float, pd.Series(x).to_list()[0].strip('][').split(', '))))
# Rename columns
floss = floss.rename(columns = lambda x : 'floss_' + str(x))
# Add back to main data frame
dfs = pd.concat([dfs[:], floss[:]], axis=1)
# Remove columns
dfs = dfs.drop(columns=['floss'])

#%%
# Take tcloss column that is a string representation of a list
# and explode to multiple columns
tcloss = dfs['tcloss'].apply(lambda x: pd.Series(map(float, pd.Series(x).to_list()[0].strip('][').split(', '))))
# Rename columns
tcloss = tcloss.rename(columns = lambda x : 'tcloss_' + str(x))
# Add back to main data frame
dfs = pd.concat([dfs[:], tcloss[:]], axis=1)
# Remove columns
dfs = dfs.drop(columns=['tcloss'])

#%%
# Take tclossc column that is a string representation of a list
# and explode to multiple columns
tclossc = dfs['tclossc'].apply(lambda x: pd.Series(map(float, pd.Series(x).to_list()[0].strip('][').split(', '))))
# Rename columns
tclossc = tclossc.rename(columns = lambda x : 'tclossc_' + str(x))
# Add back to main data frame
dfs = pd.concat([dfs[:], tclossc[:]], axis=1)
# Remove columns
dfs = dfs.drop(columns=['tclossc'])

#%%
# Group by region and scenario
dfs = dfs.groupby(['region', 'scenario'], as_index=False).sum()

#%%
# Make new columns to hold the ratio of flux lost to corridor tree cover loss
cols = ['ftcratio_{:01d}'.format(i) for i in np.arange(0, 18+1)]
# Calculate ratio
dfs[cols] = dfs.filter(regex='^floss_').div(dfs.filter(regex='^tclossc_').values)

#%%
# Melt so annual metrics plot as a line for each country
mdat = dfs.melt(['scenario', 'region'], var_name='cols',  value_name='vals')
mdat = mdat.astype({'vals': float})

# Get indices with flux loss and tree cover loss in corridor and the ratio between them 
idx = list(set([i for i in mdat.cols if any(['floss' in i, 'tclossc' in i, 'ftcratio' in i])]))
# Remove floss_0 and tclossc_0
idx = [i for i in idx if all([i != 'floss_0', i != 'tclossc_0', i != 'ftcratio_0'])]
# Get rows of flux and tree cover loss in corridors and the ratio between them
df1 = mdat[mdat['cols'].isin(idx)]
# Get neg01 scenario
#df1 = df1[df1['scenario'].isin(['neg01'])]
# Create new numeric column for flux, tree cover loss, and the ratio between them
df1 = df1.assign(yrs = np.tile(np.repeat([range(1,19,1)], np.sum(df1.cols=='floss_1')),3))

# Create label column for flux and tree cover loss in corridors and the ratio between them
def f(row):
    if 'floss' in row['cols']:
        val = 'floss'
    elif 'tclossc' in row['cols']:
        val = 'tclossc'
    elif 'ftcratio' in row['cols']:
        val = 'ftcratio'
    else:
        val = 'notype'
    return val
df1['losstype'] = df1.apply(f, axis=1)

#palette = dict(zip(dots.coherence.unique(),
#                   sns.color_palette("rocket_r", 6)))
#g = sns.relplot(x='yrs', y='vals', hue='country', data=df1, kind='line', col='region')
sns.set(style="ticks", rc={"lines.linewidth": 4, "legend.fontsize":20})
#sns.plotting_context(rc={"legend.fontsize":20})
g = sns.relplot(x='yrs', y='vals', hue='scenario', hue_order=['pos2','neg01','neg2'],
                palette=[sns.color_palette()[i] for i in [2,0,1]],
                data=df1, kind='line', col='losstype',
                row='region', col_order=['tclossc','floss','ftcratio'],
                facet_kws={'sharey': False, 'sharex': False, 'legend_out': True})

for m in range(0,3,1):
    for n in range(0,3,1):
        g.axes[m,n].tick_params(labelsize=20)
        g.axes[m,n].tick_params(labelsize=20)
        g.axes[m,n].tick_params(labelsize=20)

#g.fig.suptitle("My Title", x=0.4, y=0.98)

g.axes[0,0].set_title('', size=20)
g.axes[0,1].set_title('Africa', size=20, fontweight='bold')
g.axes[0,2].set_title('', size=20)

g.axes[1,0].set_title('', size=20)
g.axes[1,1].set_title('Asia', size=20, fontweight='bold')
g.axes[1,2].set_title('', size=20)

g.axes[2,0].set_title('', size=20)
g.axes[2,1].set_title('South America', size=20, fontweight='bold')
g.axes[2,2].set_title('', size=20)

g.fig.subplots_adjust(left=0.1, bottom=0.1, top=0.95, wspace=.4, hspace=.4)

g.axes[0,0].ticklabel_format(style='plain', axis='y')
g.axes[1,0].ticklabel_format(style='plain', axis='y')
g.axes[2,0].ticklabel_format(style='plain', axis='y')

g.axes[2,0].set_xlabel('', fontsize=20)
g.axes[2,1].set_xlabel('Loss Year', fontsize=20)
g.axes[2,2].set_xlabel('', fontsize=20)

g.axes[0,0].set_ylabel('')
g.axes[1,0].set_ylabel('Corridor Tree Cover Loss (sq km)', fontsize=20)
g.axes[2,0].set_ylabel('')

g.axes[0,1].set_ylabel('')
g.axes[1,1].set_ylabel('Corridor Flux Loss', fontsize=20)
g.axes[2,1].set_ylabel('')

g.axes[0,2].set_ylabel('')
g.axes[1,2].set_ylabel('Flux Loss/Tree Cover Loss', fontsize=20)
g.axes[2,2].set_ylabel('')

g.fig.set_size_inches(16,11)

leg = g._legend
leg.set_title('Dispersal Scenario')
new_labels = ['', 'High', 'Medium', 'Low']
for t, l in zip(leg.texts, new_labels): t.set_text(l)
    
leg.get_title().set_fontsize(20)

#%%
#------------------------------------------------------------
# SUMMARY STATS
#------------------------------------------------------------
# Sum tree cover loss in corridors by region and scenario across 1-18 indices.
# Including 0 gets the corridor pixels that didn't experience loss. Summing
# all these together actually gives total tree cover in corridors in 2000.
# Excluding the 0 index gives tree cover loss.
sdf = pd.concat([dfs.filter(items=['region','scenario']), dfs.filter(regex='^tclossc_')],
                   sort=False, axis=1)
# Drop 0th index
sdf = sdf.drop(columns=['tclossc_0'])

# Concat region and scenario to sum across years for pan-tropical tree cover loss in corridors
sdftcc = pd.concat([sdf.filter(items=['region','scenario']),sdf.sum(axis = 1, skipna = True)],
           sort=False, axis=1)
sdftcc = sdftcc.rename(columns={0: "tcloss"})
sdftcc = sdftcc.groupby(['scenario']).tcloss.sum()
print(sdftcc)

#%%
# Calculate year 2000 tree cover in corridors including ifl areas
sdf = pd.concat([dfs.filter(items=['region','scenario']), dfs.filter(regex='^tcac2000')],
                   sort=False, axis=1)
sdftcc00 = sdf.groupby(['scenario']).tcac2000.sum()
print(sdftcc00)

#%%
# Calculate year 2000 tree cover in corridors excluding ifl areas
sdf = pd.concat([dfs.filter(items=['region','scenario']), dfs.filter(regex='^tcacnoifl')],
                   sort=False, axis=1)
sdftccnoifl = sdf.groupby(['scenario']).tcacnoifl.sum()
print(sdftccnoifl)

#%%
# Test trends
import pymannkendall as mk

scen = 'neg01'
reg = 'sa'
lt = 'ftcratio'

mkt = df1[df1['scenario'].isin([scen])]
mkt = mkt[mkt['region'].isin([reg])]
mkt = mkt[mkt['losstype'].isin([lt])]
print mk.original_test(mkt.vals)
print mk.hamed_rao_modification_test(mkt.vals)
print mk.yue_wang_modification_test(mkt.vals)
print mk.sens_slope(mkt.vals)

#%%
#------------------------------------------------------------
# BY ECOREGION
#------------------------------------------------------------

#%%
# Input dir for csvs
idir = 'G:/My Drive/Projects/ifl_corridors/statsfiles/wwf_lybytc2000_fluxes/'

#%%
flist = glob.glob(os.path.join(idir, '*wwf_lybytc2000_flux_sums.csv'))
flist = [pd.read_csv(f) for f in flist]
# Concatenate files
dfs = pd.concat(flist)
# Delete extraneous columns
dfs = dfs.drop(columns=['Unnamed: 0', 'Unnamed: 0.1'])
# Rename ucode to country
dfs = dfs.rename(columns={"ucode": "ecoid"})
# Remove rows with no flux in 2000
dfs = dfs[dfs.flux2000 > 0]

#%%
# Take floss column that is a string representation of a list
# and explode to multiple columns
floss = dfs['floss'].apply(lambda x: pd.Series(map(float, pd.Series(x).to_list()[0].strip('][').replace('masked','0').split(', '))))
# Rename columns
floss = floss.rename(columns = lambda x : 'floss_' + str(x))
# Add back to main data frame
dfs = pd.concat([dfs[:], floss[:]], axis=1)
# Remove columns
dfs = dfs.drop(columns=['floss'])

#%%
# Take tcloss column that is a string representation of a list
# and explode to multiple columns
tcloss = dfs['tcloss'].apply(lambda x: pd.Series(map(float, pd.Series(x).to_list()[0].strip('][').replace('masked','0').split(', '))))
# Rename columns
tcloss = tcloss.rename(columns = lambda x : 'tcloss_' + str(x))
# Add back to main data frame
dfs = pd.concat([dfs[:], tcloss[:]], axis=1)
# Remove columns
dfs = dfs.drop(columns=['tcloss'])

#%%
# Take tclossc column that is a string representation of a list
# and explode to multiple columns
tclossc = dfs['tclossc'].apply(lambda x: pd.Series(map(float, pd.Series(x).to_list()[0].strip('][').replace('masked','0').split(', '))))
# Rename columns
tclossc = tclossc.rename(columns = lambda x : 'tclossc_' + str(x))
# Add back to main data frame
dfs = pd.concat([dfs[:], tclossc[:]], axis=1)
# Remove columns
dfs = dfs.drop(columns=['tclossc'])

#%%
# Make new columns to hold the ratio of flux lost to corridor tree cover loss
cols = ['ftcratio_{:01d}'.format(i) for i in np.arange(0, 18+1)]
# Calcualte ratio
dfs[cols] = dfs.filter(regex='^floss_').div(dfs.filter(regex='^tclossc_').values)

#%%
# Melt so annual metrics plot as a line for each country
mdat = dfs.melt(['ecoid','scenario', 'region', 'year'], var_name='cols',  value_name='vals')
mdat = mdat.astype({'vals': float})

# Get indices with flux loss and tree cover loss in corridor and the ratio between them 
idx = list(set([i for i in mdat.cols if any(['floss' in i, 'tclossc' in i, 'ftcratio' in i])]))
# Remove floss_0 and tclossc_0
idx = [i for i in idx if all([i != 'floss_0', i != 'tclossc_0', i != 'ftcratio_0'])]
# Get rows of flux and tree cover loss in corridors and the ratio between them
df1 = mdat[mdat['cols'].isin(idx)]
# Get neg01 scenario
df1 = df1[df1['scenario'].isin(['neg01'])]
# Create new numeric column for flux, tree cover loss, and the ratio between them
df1 = df1.assign(yrs = np.tile(np.repeat([range(1,19,1)], np.sum(df1.cols=='floss_1')),3))

# Create label column for flux and tree cover loss in corridors and the ratio between them
def f(row):
    if 'floss' in row['cols']:
        val = 'floss'
    elif 'tclossc' in row['cols']:
        val = 'tclossc'
    elif 'ftcratio' in row['cols']:
        val = 'ftcratio'
    else:
        val = 'notype'
    return val
df1['losstype'] = df1.apply(f, axis=1)

%matplotlib TkAgg
#palette = dict(zip(dots.coherence.unique(),
#                   sns.color_palette("rocket_r", 6)))
#g = sns.relplot(x='yrs', y='vals', hue='country', data=df1, kind='line', col='region')

g = sns.relplot(x='yrs', y='vals', hue='ecoid', data=df1, kind='line', col='losstype', row='region', facet_kws={'sharey': False, 'sharex': False})

#%%
# Test trends
import pymannkendall as mk
import collections

#scen = 'neg01'
#reg = 'sa'
#lt = ['floss', 'ftcratio', 'tclossc']

pdlist = []
for scen in ['neg01']:
    for reg in ['sa','af','as']:
        for lt in ['floss']:
            mkt = df1[df1['scenario'].isin([scen])]
            mkt = mkt[mkt['region'].isin([reg])]
            mkt = mkt[mkt['losstype'].isin([lt])]         
            ywtest = []
            sslope = []
            for t in set(mkt.ecoid.to_list()):
                try:
                    ywtest.append(mk.yue_wang_modification_test(mkt[mkt['ecoid'].isin([t])].vals).trend)
                    sslope.append(mk.sens_slope(mkt[mkt['ecoid'].isin([t])].vals))
                except ZeroDivisionError:
                    print(t)
                    print("division by zero!")
                    ywtest.append(np.nan)
                    sslope.append(np.nan)
            ywtt = collections.Counter(ywtest) 
            print(ywtt)
            print('percent of ecoregions showing increasing impact is')
            print(np.float(ywtt['increasing'])/np.sum([ywtt['no trend'],ywtt['increasing'],ywtt['decreasing']])*100.0)
            print('percent of ecoregions showing decreasing impact is')
            print(np.float(ywtt['decreasing'])/np.sum([ywtt['no trend'],ywtt['increasing'],ywtt['decreasing']])*100.0)
            print('percent of ecoregions showing no trend in impact is')
            print(np.float(ywtt['no trend'])/np.sum([ywtt['no trend'],ywtt['increasing'],ywtt['decreasing']])*100.0)
            mkt = pd.DataFrame({'ecoid': list(set(mkt.ecoid.to_list())), 'ftctrend': ywtest})
            pdlist.append(mkt)

pdlist = pd.concat(pdlist)

# Other mk options
#sns.distplot([x for x in sslope if str(x) != 'nan'])
#print mk.original_test(mkt.vals)
#print mk.hamed_rao_modification_test(mkt.vals)
#print mk.yue_wang_modification_test(mkt.vals)
#print mk.sens_slope(mkt.vals)

#%%
# Read in Ecoregion shapefile, join, and calculate trend by restoration stats
import geopandas as gpd
# Ecoregion shapes (nori indicates rock and ice ecoregion excluded)
er = gpd.read_file('G:/My Drive/Projects/common_datasets/wwf_ecoregions/Ecoregions2017/Ecoregions2017_nori.shp')

# Merge trends with ecoregions
er2 = pd.merge(er, pdlist, how='inner', on=None, left_on=['ECO_ID'], right_on=['ecoid'])
# Write to file
er2.to_file('G:/My Drive/Projects/common_datasets/wwf_ecoregions/Ecoregions2017/Ecoregions2017_nori_ftctrend.shp')

# Drop geometry column
er2 = er2.drop(columns=['geometry'])

ct = pd.crosstab(er2.NNH_NAME, er2.ftctrend)

# libraries
import matplotlib.pyplot as plt
from matplotlib import rc

# y-axis in bold
rc('font', weight='bold')
 
# Values of each group
bars1 = ct.iloc[0]
bars2 = ct.iloc[1]
bars3 = ct.iloc[2]
bars4 = ct.iloc[3]

# Heights of bars1 + bars2
barsl1 = np.add(bars1, bars2).tolist()
barsl2 = np.add(barsl1, bars3).tolist()
 
# The position of the bars on the x-axis
r = [0,1,2]
 
# Names of group and bar width
names = ['Decelerating','Accelerating','No Trend']
barWidth = 0.75
 
# Create brown bars
p1 = plt.bar(r, bars1, color='#7f6d5f', edgecolor='white', width=barWidth)
# Create green bars (middle), on top of the firs ones
p2 = plt.bar(r, bars2, bottom=bars1, color='#557f2d', edgecolor='white', width=barWidth)
# Create green bars (top)
p3 = plt.bar(r, bars3, bottom=barsl1, color='#2d7f5e', edgecolor='white', width=barWidth)
# Create green bars (top)
p4 = plt.bar(r, bars4, bottom=barsl2, color='cornflowerblue', edgecolor='white', width=barWidth)
     
# Custom X axis
plt.xticks(r, names, fontweight='bold', fontsize=20)
plt.xlabel("Trend in Deforestation's Impact on Connectivity", fontsize=20, weight='bold')
plt.yticks(fontsize=20)

plt.legend((p4[0], p3[0], p2[0], p1[0]), ('NI', 'NCR', 'NCRHP', 'HP'), fontsize=20)

plt.title('Ecoregion Status and Connectivity Trends', fontsize=20, weight='bold')
# Show graphic
plt.show()


#%%
#------------------------------------------------------------
# BY IFL ZONE
#------------------------------------------------------------

#%%
# Input dir for csvs
idir = 'G:/My Drive/Projects/ifl_corridors/statsfiles/iflzone_lybytc2000_fluxes/'

#%%
flist = glob.glob(os.path.join(idir, '*_iflzone_lybytc2000_flux_sums.csv'))
flist = [pd.read_csv(f) for f in flist]
# Concatenate files
dfs = pd.concat(flist)
# Delete extraneous columns
dfs = dfs.drop(columns=['Unnamed: 0', 'Unnamed: 0.1'])
# Rename ucode to country
dfs = dfs.rename(columns={"ucode": "iflid"})
# Remove rows with no flux in 2000
dfs = dfs[dfs.flux2000 > 0]

#%%
# Take floss column that is a string representation of a list
# and explode to multiple columns
floss = dfs['floss'].apply(lambda x: pd.Series(map(float, pd.Series(x).to_list()[0].strip('][').replace('masked','0').split(', '))))
# Rename columns
floss = floss.rename(columns = lambda x : 'floss_' + str(x))
# Add back to main data frame
dfs = pd.concat([dfs[:], floss[:]], axis=1)
# Remove columns
dfs = dfs.drop(columns=['floss'])

#%%
# Take tcloss column that is a string representation of a list
# and explode to multiple columns
tcloss = dfs['tcloss'].apply(lambda x: pd.Series(map(float, pd.Series(x).to_list()[0].strip('][').replace('masked','0').split(', '))))
# Rename columns
tcloss = tcloss.rename(columns = lambda x : 'tcloss_' + str(x))
# Add back to main data frame
dfs = pd.concat([dfs[:], tcloss[:]], axis=1)
# Remove columns
dfs = dfs.drop(columns=['tcloss'])

#%%
# Take tclossc column that is a string representation of a list
# and explode to multiple columns
tclossc = dfs['tclossc'].apply(lambda x: pd.Series(map(float, pd.Series(x).to_list()[0].strip('][').replace('masked','0').split(', '))))
# Rename columns
tclossc = tclossc.rename(columns = lambda x : 'tclossc_' + str(x))
# Add back to main data frame
dfs = pd.concat([dfs[:], tclossc[:]], axis=1)
# Remove columns
dfs = dfs.drop(columns=['tclossc'])

#%%
# Make new columns to hold the ratio of flux lost to corridor tree cover loss
cols = ['ftcratio_{:01d}'.format(i) for i in np.arange(0, 18+1)]
# Calcualte ratio
dfs[cols] = dfs.filter(regex='^floss_').div(dfs.filter(regex='^tclossc_').values)

#%%
# Melt so annual metrics plot as a line for each country
mdat = dfs.melt(['iflid','scenario', 'region', 'year'], var_name='cols',  value_name='vals')
mdat = mdat.astype({'vals': float})

# Get indices with flux loss and tree cover loss in corridor and the ratio between them 
idx = list(set([i for i in mdat.cols if any(['floss' in i, 'tclossc' in i, 'ftcratio' in i])]))
# Remove floss_0 and tclossc_0
idx = [i for i in idx if all([i != 'floss_0', i != 'tclossc_0', i != 'ftcratio_0'])]
# Get rows of flux and tree cover loss in corridors and the ratio between them
df1 = mdat[mdat['cols'].isin(idx)]
# Create new numeric column for flux, tree cover loss, and the ratio between them
df1 = df1.assign(yrs = np.tile(np.repeat([range(1,19,1)], np.sum(df1.cols=='floss_1')),3))

# Create label column for flux and tree cover loss in corridors and the ratio between them
def f(row):
    if 'floss' in row['cols']:
        val = 'floss'
    elif 'tclossc' in row['cols']:
        val = 'tclossc'
    elif 'ftcratio' in row['cols']:
        val = 'ftcratio'
    else:
        val = 'notype'
    return val
df1['losstype'] = df1.apply(f, axis=1)

#----
# Plot trends for individual IFLs
#sa301 = df1[(df1.iflid == 301) & (df1.region == 'sa')]
sa153 = df1[(df1.iflid == 153) & (df1.region == 'sa')]

sns.set(style="ticks", rc={"lines.linewidth": 4, "legend.fontsize":20})

g = sns.relplot(x='yrs', y='vals', hue='scenario', data=sa153, kind='line', col='losstype', facet_kws={'sharey': False, 'sharex': False})

for n in range(0,3,1):
    g.axes[0,n].tick_params(labelsize=20)
    g.axes[0,n].tick_params(labelsize=20)
    g.axes[0,n].tick_params(labelsize=20)

g.axes[0,0].set_title('', size=20)
g.axes[0,1].set_title('', size=20, fontweight='bold')
g.axes[0,2].set_title('', size=20)

g.fig.subplots_adjust(left=0.1, bottom=0.2, top=0.95, wspace=.4, hspace=.4)

g.axes[0,0].ticklabel_format(style='plain', axis='y')
#g.axes[1,0].ticklabel_format(style='plain', axis='y')
#g.axes[2,0].ticklabel_format(style='plain', axis='y')

g.axes[0,0].set_xlabel('', fontsize=20)
g.axes[0,1].set_xlabel('Loss Year', fontsize=20)
g.axes[0,2].set_xlabel('', fontsize=20)

g.axes[0,0].set_ylabel('Corridor Tree Cover Loss (sq km)', fontsize=20)
#g.axes[1,0].set_ylabel('Corridor Tree Cover Loss (sq km)', fontsize=20)
#g.axes[2,0].set_ylabel('')

g.axes[0,1].set_ylabel('Corridor Flux Loss', fontsize=20)
#g.axes[1,1].set_ylabel('Corridor Flux Loss', fontsize=20)
#g.axes[2,1].set_ylabel('')

g.axes[0,2].set_ylabel('Flux Loss/Tree Cover Loss', fontsize=20)
#g.axes[1,2].set_ylabel('Flux Loss/Tree Cover Loss', fontsize=20)
#g.axes[2,2].set_ylabel('')

g.fig.set_size_inches(20,6)

leg = g._legend
leg.set_title('Dispersal Scenario')
new_labels = ['', 'High', 'Medium', 'Low']
for t, l in zip(leg.texts, new_labels): t.set_text(l)
    
leg.get_title().set_fontsize(20)

#%%
# Test trends
import pymannkendall as mk
import collections

pdlist = []
for scen in ['neg2','neg01','pos2']:
    for reg in ['sa','af','as']:
        for lt in ['floss', 'ftcratio', 'tclossc']:
            mkt = df1[df1['scenario'].isin([scen])]
            mkt = mkt[mkt['region'].isin([reg])]
            mkt = mkt[mkt['losstype'].isin([lt])]         
            ywtest = []
            sslope = []
            for t in set(mkt.iflid.to_list()):
                try:
                    ywtest.append(mk.yue_wang_modification_test(mkt[mkt['iflid'].isin([t])].vals).trend)
                    sslope.append(mk.sens_slope(mkt[mkt['iflid'].isin([t])].vals))
                except ZeroDivisionError:
                    print(t)
                    print("division by zero!")
                    ywtest.append(np.nan)
                    sslope.append(np.nan)
            ywtt = collections.Counter(ywtest) 
            print(ywtt)
            print('percent of zones showing increasing impact is')
            print(np.float(ywtt['increasing'])/np.sum([ywtt['no trend'],ywtt['increasing'],ywtt['decreasing']])*100.0)
            print('percent of zones showing decreasing impact is')
            print(np.float(ywtt['decreasing'])/np.sum([ywtt['no trend'],ywtt['increasing'],ywtt['decreasing']])*100.0)
            print('percent of zones showing no trend in impact is')
            print(np.float(ywtt['no trend'])/np.sum([ywtt['no trend'],ywtt['increasing'],ywtt['decreasing']])*100.0)
            mkt = pd.DataFrame({'iflid': list(set(mkt.iflid.to_list())), 'trendir': ywtest,
                                'trendsl': sslope, 'losstype': lt, 'region': reg, 'scenario': scen})
            pdlist.append(mkt)

pdlist = pd.concat(pdlist)

# Get top X% of flux loss increases
for r in ['sa', 'af', 'as']:
    for s in ['neg2', 'neg01', 'pos2']:
        adat = pdlist[(pdlist['region']==r) & (pdlist['scenario']==s) & (pdlist['losstype']=='floss')]
        topq = np.quantile(adat.trendsl[adat.trendsl >=0], 0.95)
        topqs = set(adat[adat.trendsl >= topq].iflid.tolist())
        adat = pdlist[(pdlist['region']==r) & (pdlist['scenario']==s) & (pdlist['losstype']=='ftcratio')]
        topq = np.quantile(adat.trendsl[adat.trendsl >=0], 0.95)
        topqs2 = set(adat[adat.trendsl >= topq].iflid.tolist())
        print(topqs & topqs2)

#----
# Create dataframe for each set of trends
for r in ['sa', 'af', 'as']:
    for s in ['neg2', 'neg01', 'pos2']:
        d1 = pd.merge(pdlist[(pdlist['region']==r) & (pdlist['scenario']==s) & (pdlist['losstype']=='floss')],
                        pdlist[(pdlist['region']==r) & (pdlist['scenario']==s) & (pdlist['losstype']=='ftcratio')],
                        left_on='iflid', right_on='iflid')
        d2 = pd.merge(d1, pdlist[(pdlist['region']==r) & (pdlist['scenario']==s) & (pdlist['losstype']=='tclossc')],
                        left_on='iflid', right_on='iflid')
        # Drop columns
        d2 = d2.drop(columns=['losstype_x', 'losstype_y', 'region_y', 'scenario_y', 'losstype', 'region', 'scenario'])
        # Rename columns
        d2 = d2.rename(columns={'region_x': 'region', 'scenario_x': 'scenario', 'trendir_x': 'tdfloss',
                       'trendsl_x': 'tsfloss', 'trendir_y': 'tdftc', 'trendsl_y': 'tsftc',
                       'trendir': 'tdtcloss', 'trendsl': 'tstcloss'})
        d2['tgroups'] = d2.tdfloss + '_' + d2.tdftc + '_' + d2.tdtcloss
        d2['fgroups'] = d2.tdfloss + '_' + d2.tdftc
        #print(collections.Counter(d2.tgroups))
        # Write to file
        d2.to_csv(idir + 'df' + s + '_' + r + '.csv')


#----
# Barplots of trends by ifl zone
bplist = []
for r in ['sa', 'af', 'as']:
    for s in ['neg01']:
        d1 = pd.merge(pdlist[(pdlist['region']==r) & (pdlist['scenario']==s) & (pdlist['losstype']=='floss')],
                        pdlist[(pdlist['region']==r) & (pdlist['scenario']==s) & (pdlist['losstype']=='ftcratio')],
                        left_on='iflid', right_on='iflid')
        d2 = pd.merge(d1, pdlist[(pdlist['region']==r) & (pdlist['scenario']==s) & (pdlist['losstype']=='tclossc')],
                        left_on='iflid', right_on='iflid')
        # Drop columns
        d2 = d2.drop(columns=['losstype_x', 'losstype_y', 'region_y', 'scenario_y', 'losstype', 'region', 'scenario'])
        # Rename columns
        d2 = d2.rename(columns={'region_x': 'region', 'scenario_x': 'scenario', 'trendir_x': 'tdfloss',
                       'trendsl_x': 'tsfloss', 'trendir_y': 'tdftc', 'trendsl_y': 'tsftc',
                       'trendir': 'tdtcloss', 'trendsl': 'tstcloss'})
        d2['tgroups'] = d2.tdfloss + '_' + d2.tdftc + '_' + d2.tdtcloss
        d2['fgroups'] = d2.tdfloss + '_' + d2.tdftc
        d2.fgroups = d2.fgroups.str.replace('no trend', 'steady', regex=True)
        
        ct = pd.crosstab(d2.iflid, d2.fgroups)
        ct = ct.sum(axis = 0, skipna = True)
        corder = ['decreasing_decreasing', 'decreasing_steady', 'decreasing_increasing',
          'steady_decreasing', 'steady_steady', 'steady_increasing',
          'increasing_decreasing', 'increasing_steady', 'increasing_increasing']
        ct = ct[corder]
        ct = pd.DataFrame({'trend': list(ct.index), 'count': list(ct), 'region': r, 'scenario': s})
        ct = ct.fillna(0)
        ct = ct.astype({'count': 'int16'})
        bplist.append(ct)

# Make dataframe from trend count dataframes
bpdf = pd.DataFrame({'trend': corder, 'sa': bplist[0]['count'], 'af': bplist[1]['count'], 'as': bplist[2]['count']})
# Reorder columns
bpdf = bpdf[['trend', 'sa', 'af', 'as']]
bpdf.sum(skipna=True)
pctrend = bpdf.loc[:,['sa', 'af', 'as']].div(bpdf.sum(axis=0, skipna=True))
pctrend = pctrend*100

# Increasing trend in flux loss
print(np.sum(pctrend['sa'][6:9]))
print(np.sum(pctrend['af'][6:9]))
print(np.sum(pctrend['as'][6:9]))

# Increasing trend in connectivity impact
print(np.sum(pctrend['sa'][[2,5,8]]))
print(np.sum(pctrend['af'][[2,5,8]]))
print(np.sum(pctrend['as'][[2,5,8]]))

# libraries
import matplotlib.pyplot as plt
from matplotlib import rc

# y-axis in bold
rc('font', weight='bold')
 
bars1 = bpdf.loc[8, ['sa','af','as']]
bars2 = bpdf.loc[7, ['sa','af','as']]
bars3 = bpdf.loc[6, ['sa','af','as']]
bars4 = bpdf.loc[5, ['sa','af','as']]
bars5 = bpdf.loc[4, ['sa','af','as']]
bars6 = bpdf.loc[3, ['sa','af','as']]
bars7 = bpdf.loc[2, ['sa','af','as']]
bars8 = bpdf.loc[1, ['sa','af','as']]
bars9 = bpdf.loc[0, ['sa','af','as']]

# Heights of bars1 + bars2
barsl1 = np.add(bars1, bars2).tolist()
barsl2 = np.add(barsl1, bars3).tolist()
barsl3 = np.add(barsl2, bars4).tolist()
barsl4 = np.add(barsl3, bars5).tolist()
barsl5 = np.add(barsl4, bars6).tolist()
barsl6 = np.add(barsl5, bars7).tolist()
barsl7 = np.add(barsl6, bars8).tolist()
barsl8 = np.add(barsl7, bars9).tolist()
 
# The position of the bars on the x-axis
r = [0,1,2]
r = [0,.75,1.5]
 
# Names of group and bar width
names = ['SA','AF','AS']
barWidth = 0.5

cpal = sns.color_palette("Oranges", n_colors=9)
cpal2 = sns.color_palette("Reds", n_colors=9)

# Stack bars

p2 = plt.bar(r, bars2, bottom=bars1, color=cpal2[7], edgecolor='white', width=barWidth)
p3 = plt.bar(r, bars3, bottom=barsl1, color=cpal2[6], edgecolor='white', width=barWidth)

p5 = plt.bar(r, bars5, bottom=barsl3, color=cpal[4], edgecolor='white', width=barWidth)
p6 = plt.bar(r, bars6, bottom=barsl4, color=cpal[3], edgecolor='white', width=barWidth)

p8 = plt.bar(r, bars8, bottom=barsl6, color=cpal[1], edgecolor='white', width=barWidth)
p9 = plt.bar(r, bars9, bottom=barsl7, color=cpal[0], edgecolor='white', width=barWidth)     

p1 = plt.bar(r, bars1, color=cpal2[8], edgecolor='black', width=barWidth)
p4 = plt.bar(r, bars4, bottom=barsl2, color=cpal[5], edgecolor='black', width=barWidth)
p7 = plt.bar(r, bars7, bottom=barsl5, color=cpal[2], edgecolor='black', width=barWidth)

# Custom X axis
plt.xticks(r, names, fontweight='bold', fontsize=20)
plt.xlabel("Region", fontsize=20, weight='bold')
plt.yticks(fontsize=20)

plt.legend((p9[0], p8[0], p7[0], p6[0], p5[0], p4[0], p3[0], p2[0], p1[0]),
           ('DD', 'DS', 'DI', 'SD', 'SS', 'SI', 'ID', 'IS', 'II'), fontsize=20)

plt.title('Trends in Connectivity Loss and Connectivity Impact of Deforestation by IFL', fontsize=20, weight='bold')
# Show graphic
plt.show()


#----
# Compare trends with tree cover stats in each IFL
#%%
# Input dir for csvs
idir = 'G:/My Drive/Projects/ifl_corridors/statsfiles/iflzone_lybytc2000_fluxes/'

#%%
flist = glob.glob(os.path.join(idir, '*_iflzone_lybytc2000_flux_sums.csv'))
flist = [pd.read_csv(f) for f in flist]
# Concatenate files
dfs = pd.concat(flist)
# Delete extraneous columns
dfs = dfs.drop(columns=['Unnamed: 0', 'Unnamed: 0.1'])
# Rename ucode to country
dfs = dfs.rename(columns={"ucode": "iflid"})
# Remove rows with no flux in 2000
dfs = dfs[dfs.flux2000 > 0]
# Drop columns
dfs = dfs.drop(columns=['floss', 'tcloss', 'tclossc', 'year'])

trendtclist = []
for r in ['sa', 'af', 'as']:
    for s in ['neg01']:
        d1 = pd.merge(pdlist[(pdlist['region']==r) & (pdlist['scenario']==s) & (pdlist['losstype']=='floss')],
                        pdlist[(pdlist['region']==r) & (pdlist['scenario']==s) & (pdlist['losstype']=='ftcratio')],
                        left_on='iflid', right_on='iflid')
        d2 = pd.merge(d1, pdlist[(pdlist['region']==r) & (pdlist['scenario']==s) & (pdlist['losstype']=='tclossc')],
                        left_on='iflid', right_on='iflid')
        # Drop columns
        d2 = d2.drop(columns=['losstype_x', 'losstype_y', 'region_y', 'scenario_y', 'losstype', 'region', 'scenario'])
        # Rename columns
        d2 = d2.rename(columns={'region_x': 'region', 'scenario_x': 'scenario', 'trendir_x': 'tdfloss',
                       'trendsl_x': 'tsfloss', 'trendir_y': 'tdftc', 'trendsl_y': 'tsftc',
                       'trendir': 'tdtcloss', 'trendsl': 'tstcloss'})
        d2['tgroups'] = d2.tdfloss + '_' + d2.tdftc + '_' + d2.tdtcloss
        d2['fgroups'] = d2.tdfloss + '_' + d2.tdftc
        d2.fgroups = d2.fgroups.str.replace('no trend', 'steady', regex=True)
        # Merge with data on tree and land area
        d3 = pd.merge(d2, dfs[(dfs['region']==r) & (dfs['scenario']==s)], left_on='iflid', right_on='iflid')
        # Calculate percent tree cover in zone and in corridors in zone
        d3['tcap2000'] = d3.tca2000/d3.landarea*100
        d3['tcacp2000'] = d3.tcac2000/d3.landareac*100
        # Drop region_y and scenario_y columns
        d3 = d3.drop(columns=['region_y', 'scenario_y'])
        # Rename columns
        d3 = d3.rename(columns={'region_x': 'region', 'scenario_x': 'scenario'})
        # Append to list
        trendtclist.append(d3)

tdf = pd.concat(trendtclist)
# Add column for square root of flux 2000
tdf['flux2000sqrt'] = np.sqrt(tdf['flux2000'])

# Plot with regions as columns
#cpal2 = sns.color_palette("Reds", n_colors=9)
# Create palette and reverse
cpal = sns.color_palette("Oranges", n_colors=9)
cpal.reverse()
# Order of group plotting
ho = ['increasing_increasing','increasing_steady','steady_increasing','increasing_decreasing','decreasing_increasing','steady_steady','decreasing_steady','steady_decreasing','decreasing_decreasing']
g = sns.relplot(x='tsfloss', y='tsftc', hue='fgroups', size='flux2000sqrt',
                sizes=(10, 2000), data=tdf, kind='scatter', col='region',
                row='scenario', facet_kws={'sharey': False, 'sharex': False},
                hue_order=ho, palette=cpal)

# Add zero lines
for i in range(0,3):
    g.axes[0,i].axhline(0,ls='--', color='black')
    g.axes[0,i].axvline(0,ls='--', color='black')
    g.axes[0,i].tick_params(labelsize=20)
# Set panel top labels
g.axes[0,0].set_title('South America', size=20, fontweight='bold')
g.axes[0,1].set_title('Africa', size=20, fontweight='bold')
g.axes[0,2].set_title('Asia', size=20, fontweight='bold')
# Set panel x labels
g.axes[0,0].set_xlabel('')
g.axes[0,1].set_xlabel('Connectivity Loss', size=20, fontweight='bold')
g.axes[0,2].set_xlabel('')
# Set panel y labels
g.axes[0,0].set_ylabel('Connectivity Marginal Impact', size=20, fontweight='bold')
g.axes[0,1].set_ylabel('')
g.axes[0,2].set_ylabel('')

g._legend.texts[0].set_text('Trend Groups')
g._legend.texts[10].set_text('Year 2000 Flux')

for i,lh in enumerate(g._legend.legendHandles):
    if i < 10:
        lh._sizes = [75] 

g.figure.set_size_inches(16,11)

plt.legend((p9[0], p8[0], p7[0], p6[0], p5[0], p4[0], p3[0], p2[0], p1[0]),
           ('DD', 'DS', 'DI', 'SD', 'SS', 'SI', 'ID', 'IS', 'II'), fontsize=20)

plt.title('Trends in Connectivity Loss and Connectivity Impact of Deforestation by IFL', fontsize=20, weight='bold')
# Show graphic
plt.show()


sattc = trendtclist[0]
g = sns.scatterplot(sattc.tsfloss, sattc.tsftc, size=sattc.tcac2000, sizes=(100, 1000), hue=sattc.fgroups)
g.axes.axhline(0,ls='--', color='black')
g.axes.axvline(0,ls='--', color='black')
g.axes.set_xlabel('Trend in Connectivity Loss',fontsize=30)
g.axes.set_ylabel('Trend in Marginal Connectivity Impact',fontsize=30)
g.axes.tick_params(labelsize=30)



g = sns.relplot(x='yrs', y='vals', hue='scenario', data=sa301, kind='line', col='losstype', facet_kws={'sharey': False, 'sharex': False})




afttc = trendtclist[1]
g = sns.scatterplot(afttc.tsfloss, afttc.tsftc, size=afttc.flux2000, sizes=(100, 5000), hue=afttc.fgroups)
g.axes.axhline(0,ls='--', color='black')
g.axes.axvline(0,ls='--', color='black')
g.axes.set_xlabel('Trend in Connectivity Loss',fontsize=30)
g.axes.set_ylabel('Trend in Marginal Connectivity Impact',fontsize=30)
g.axes.tick_params(labelsize=30)

asttc = trendtclist[2]
g = sns.scatterplot(asttc.tsfloss, asttc.tsftc, size=asttc.flux2000, sizes=(100, 2000), hue=asttc.fgroups)
g.axes.axhline(0,ls='--', color='black')
g.axes.axvline(0,ls='--', color='black')
g.axes.set_xlabel('Trend in Connectivity Loss',fontsize=30)
g.axes.set_ylabel('Trend in Marginal Connectivity Impact',fontsize=30)
g.axes.tick_params(labelsize=30)



sattc = trendtclist[0]
g = sns.scatterplot(sattc.flux2000**0.5, sattc.tsfloss**0.5,  size=sattc.tsftc, sizes=(1, 1000), hue=sattc.fgroups)
g.axes.axhline(0,ls='--', color='black')
g.axes.axvline(0,ls='--', color='black')
g.axes.set_xlabel('Year 2000 Flux', fontsize=30)
g.axes.set_ylabel('Trend in Connectivity Loss', fontsize=30)
g.axes.tick_params(labelsize=30)





sns.scatterplot(atest.tca2000,atest.tsftc, hue=atest.tstcloss)
sns.scatterplot(atest.tstcloss,atest.tsftc, hue=atest.tstcloss)

sns.scatterplot(atest.tca2000,atest.tstcloss, hue=atest.tclass)

atest['tclass'] = atest.tdfloss + '_' + atest.tdftc

afttc.groupby('fgroups', as_index=False)['tcac2000'].mean()

atest.groupby('tdfloss', as_index=False)['tcap2000'].mean()

atest.groupby('tdftc', as_index=False)['tcap2000'].hist()

atest[atest['tdftc']=='no trend']['tcap2000'].hist()
atest[atest['tdftc']=='decreasing']['tcap2000'].hist()
atest[atest['tdftc']=='increasing']['tcap2000'].hist()

atest[atest['tdfloss']=='no trend']['tcap2000'].hist()
atest[atest['tdfloss']=='decreasing']['tcap2000'].hist()
atest[atest['tdfloss']=='increasing']['tcap2000'].hist()

atest.to_csv(idir + 'dfneg01_sa_trend_by_treecover.csv')
