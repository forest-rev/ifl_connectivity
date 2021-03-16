# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 20:32:37 2019
Make a bubble plot of country fluxes
@author: pj276
"""

#%%
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib as plt
import glob
from decimal import *

#%matplotlib qt

#%%
# List files
flist = glob.glob('G:/My Drive/Projects/ifl_corridors/statsfiles/gadm_fluxes/*')
flist = [pd.read_csv(f) for f in flist]
# Concatenate files
dfs = pd.concat(flist)
# Delete extraneous columns
dfs = dfs.drop(columns=['Unnamed: 0', 'Unnamed: 0.1'])
# Rename ucode to country
dfs = dfs.rename(columns={"ucode": "country"})

#%%
# Graph only neg01 scenario
gnew = dfs[dfs['scenario'] == 'neg01']

df1 = gnew[gnew['year'] == 2000]
df2 = gnew[gnew['year'] == 2018]

#%%
# Join them together
gnew = df1.set_index('country').join(df2.set_index('country'),lsuffix='_1', rsuffix='_2')
# Calculate flux percent change
gnew['fluxpc'] = np.abs((gnew['flux_2'] - gnew['flux_1'])/gnew['flux_1']*100)
# Calculate flux change
gnew['fluxch'] = (gnew['flux_2'] - gnew['flux_1'])
# Calculate flux change share
gnew['fluxchsh'] = np.abs(gnew['fluxch'])/np.sum(gnew['flux_1'])

# Remove countries with low flux in 2000
gnew = gnew[gnew['flux_1'] > 10]

# Calculate flux share per country in 2000
gnew['fluxsh2000'] = gnew['flux_1']/np.sum(gnew['flux_1'])*100

# Calculate log of flux in 2000
gnew['flux_1_ln'] = np.log(gnew['flux_1'])

# Get country codes
gnew['codes'] = gnew.index

# Set hue variable
gnew['fluxchhue'] = np.abs(gnew['flux_2'] - gnew['flux_1'])

gnew.loc[gnew['fluxch'] < 0, 'fluxchhue'] = np.log(gnew['fluxchhue'])

#%%
# Production bubble plot
#cmap = sns.light_palette("orange",as_cmap=True)
%matplotlib qt

# Get only flux loss
gnewneg = gnew.loc[gnew['fluxch'] < 0]

# Get absolute value of flux change
gnewneg.loc[:,'fluxchpos'] = np.abs(gnewneg['fluxch'])
# Get square root of flux change and transform to use as proportional to area
gnewneg.loc[:,'fchdiam'] = np.sqrt(gnewneg['fluxchpos'])

# Bubble plot
g = sns.scatterplot('fluxpc', 'flux_1', size='fchdiam', sizes=(1, 7920),
                     data=gnewneg, alpha=0.75, legend=False)
mylab = str(np.mean(gnewneg['fluxpc'][gnewneg['fluxch'] < 0]).round(2))
g.axes.text(x=np.mean(gnewneg['fluxpc'][gnewneg['fluxch'] < 0]), y=10, s=mylab, horizontalalignment='center', size='large', color='green', weight='semibold')
g.axes.set_xlabel('Reduction in IFL Flux 2000-2018 (%)',fontsize=30)
g.axes.set_ylabel('IFL Flux in 2000',fontsize=30)
g.axes.tick_params(labelsize=30)
g.axes.set_yscale('log')
g.figure.set_size_inches(16,11)
g.axes.set_ylim(bottom=None,top=10**10)

#For each point, we add a text inside the bubble
for line in range(0,gnewneg.shape[0]):
     g.axes.text(gnewneg.fluxpc[line], gnewneg.flux_1[line], gnewneg.codes[line], horizontalalignment='center', size='medium', color='black', weight='semibold')

# Add vertical line for mean of percent flux loss
g.axes.axvline(np.mean(gnewneg['fluxpc'][gnewneg['fluxch'] < 0]), color='green', linestyle='--')

# List is representative values of absolute flux change
sizelist = [900, 90000, 9000000, 81000000]
#sizelist = [1000, 30000, 300000, 3000000, 60726400]
# List is the square root of representative values of absolute flux change
for m, area in enumerate([30, 300, 3000, 9000]):
#for m, area in enumerate([30, 170, 550, 7800]):
    g.scatter([], [], c='k', alpha=0.3, s=area,
                label="{:.2E}".format(Decimal(sizelist[m])))
legend = g.legend(scatterpoints=1, frameon=False, fontsize=20,
                  labelspacing=1, title='IFL Flux Change', handletextpad=2)
legend.get_title().set_fontsize(20)
# Save figure
g.figure.savefig('G:/My Drive/Projects/ifl_corridors/figs/flux_by_country_change.png')

#%% 
# Read in ifl status from 2017 paper
ifl = pd.read_csv('G:/My Drive/Projects/ifl_corridors/statsfiles/ifl_stats.csv', encoding='cp1252')
ifl = ifl.merge(gnewneg, left_on='Country code', right_on='codes', how='inner')
ifl['ifl_absarea_red'] = ifl['IFL area red']/100*ifl['IFL 2000']

h = sns.scatterplot('fluxchpos', 'ifl_absarea_red', size='fchdiam', sizes=(1, 7920),
                    data=ifl, alpha=0.75, legend=False, color='green')
h.axes.set_xscale('log')
h.axes.set_yscale('log')
h.axes.set_xlim(left=None,right=10**8.5)
h.axes.set_ylim(bottom=None,top=10**5)

for line in range(0,ifl.shape[0]):
     h.axes.text(ifl.fluxchpos[line], ifl['ifl_absarea_red'][line], ifl.codes[line], horizontalalignment='center', size='medium', color='black', weight='semibold')

# Add vertical line for mean of percent flux loss
h.axes.axvline(np.mean(ifl['fluxchpos']), color='green', linestyle='--')
# Add horizontal line for mean of percent ifl loss
h.axes.axhline(np.mean(ifl['ifl_absarea_red']), color='green', linestyle='--')
# Increas font sizes
h.axes.set_xlabel('Reduction in IFL Flux 2000-2012',fontsize=30)
h.axes.set_ylabel('Reduction in IFL Area 2000-2013',fontsize=30)
h.axes.tick_params(labelsize=30)
h.figure.set_size_inches(16,11)
h.figure.savefig('G:/My Drive/Projects/ifl_corridors/figs/ifl_flux_comp_by_country.png')
# Get correlation between absolute ifl and flux change
np.corrcoef(ifl.fluxchpos,ifl['ifl_absarea_red']) # 0.9
np.corrcoef(np.log(ifl.fluxchpos),np.log(ifl['ifl_absarea_red'])) # 0.56
# Get correlation coefficient between ifl area in 2000 and flux in 2000
np.corrcoef(ifl['IFL 2000'],ifl['flux_1']) # 0.9
