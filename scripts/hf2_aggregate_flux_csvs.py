# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 11:29:07 2019
This script concatenates csvs from zone flux summaries and writes to file
@author: pj276
"""

#%%
import pandas as pd
import sys, glob
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
import param_file as p1

#%%
# Top level dir
idir = p1.idir #'/scratch/pj276/ifl_corridors'

# Flux directory
fluxdir = p1.fluxdir #'flux_spp'

# Output flux name
fluxcsv = p1.fluxcsv #'spp_flux_sums.csv'

# Directory holding output aggregated csvs
aggcsvdir = p1.aggcsvdir

# Suffix to filter file results
endpatt = p1.endpatt

#%%
# Region
i = p1.region #['af','as','sa']
# Scenario
j = p1.scenario #['neg2', 'neg01', 'pos2']
# Year
#year = ['2000', '2013', '2018']
year = p1.yearlist
# Folder year equals tiff year switch
tiffequalsfolder = p1.tiffequalsfolder

# Empty list to hold list of flux dataframes
dflist = []
if tiffequalsfolder == False:
    for k in year:
        if k == '2000':
            kk = '2000'
        if k == '2013':
            kk = '2013'
        if k == '2018':
            kk = '2013'
        flist = glob.glob(idir + '/rp5k130k90m_' + j + '_' + i + '_' + kk + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '*')
        flist = [m for m in flist if m.endswith(endpatt)]
        dfs = [pd.read_csv(m) for m in flist]
        dfs = pd.concat(dfs, ignore_index=True)
        #dfs = dfs[dfs['flux'] > 0]
        #dfs['fluxf'] = dfs['flux']/np.sum(dfs['flux'])*100
        dflist.append(dfs)
    dfs = pd.concat(dflist, ignore_index=True)
    dfs.to_csv(idir + '/' + aggcsvdir + '/' + i + '_' + j + '_' + k + '_' + fluxcsv)

if tiffequalsfolder == True:
    for k in year:
        if k == '2000':
            kk = '2000'
        if k == '2013':
            kk = '2013'
        if k == '2018':
            kk = '2018'
        flist = glob.glob(idir + '/rp5k130k90m_' + j + '_' + i + '_' + kk + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '*')
        flist = [m for m in flist if m.endswith(endpatt)]
        dfs = [pd.read_csv(m) for m in flist]
        dfs = pd.concat(dfs, ignore_index=True)
        #dfs = dfs[dfs['flux'] > 0]
        #dfs['fluxf'] = dfs['flux']/np.sum(dfs['flux'])*100
        dflist.append(dfs)
    dfs = pd.concat(dflist, ignore_index=True)
    dfs.to_csv(idir + '/' + aggcsvdir + '/' + i + '_' + j + '_' + k + '_' + fluxcsv)