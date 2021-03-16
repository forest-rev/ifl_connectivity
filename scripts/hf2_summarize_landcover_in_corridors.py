# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 08:23:12 2019
Link up ESA 300m land cover names with codes from raster attribute tables.
@author: pj276
"""
#---------------------------------
import pandas as pd, numpy as np
import seaborn as sns

#---------------------------------
# Raster attribute table names
# Africa
aflc = "G:/My Drive/Projects/ifl_corridors/statsfiles/lc_stats/esa_2012_lc_in_af_neg01_corridors_no_ifl.csv"
# Asia
aslc = "G:/My Drive/Projects/ifl_corridors/statsfiles/lc_stats/esa_2012_lc_in_as_neg01_corridors_no_ifl.csv"
# South America
salc = "G:/My Drive/Projects/ifl_corridors/statsfiles/lc_stats/esa_2012_lc_in_sa_neg01_corridors_no_ifl.csv"

#---------------------------------
# Set up land cover class names and summarize land cover in corridor matrix by region
lcodes = [10,11,12,20,30,40,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,150,151,152,153,160,170,180,190,200,201,202,210,220]
lnames = ['Cropland, rainfed','Herbaceous cover','Tree or shrub cover','Cropland, irrigated or post-flooding',
          'Mosaic cropland (>50%) / natural vegetation (tree, shrub, herbaceous cover) (<50%)','Mosaic natural vegetation (tree, shrub, herbaceous cover) (>50%) / cropland (<50%)',
          'Tree cover, broadleaved, evergreen, closed to open (>15%)','Tree cover, broadleaved, deciduous, closed to open (>15%)','Tree cover, broadleaved, deciduous, closed (>40%)',
          'Tree cover, broadleaved, deciduous, open (15‐40%)','Tree cover, needleleaved, evergreen, closed to open (>15%)','Tree cover, needleleaved, evergreen, closed (>40%)',
          'Tree cover, needleleaved, evergreen, open (15‐40%)','Tree cover, needleleaved, deciduous, closed to open (>15%)',
          'Tree cover, needleleaved, deciduous, closed (>40%)','Tree cover, needleleaved, deciduous, open (15‐40%)',
          'Tree cover, mixed leaf type (broadleaved and needleleaved)','Mosaic tree and shrub (>50%) / herbaceous cover (<50%)',
          'Mosaic herbaceous cover (>50%) / tree and shrub (<50%)','Shrubland','Evergreen shrubland','Deciduous shrubland',
          'Grassland','Lichens and mosses','Sparse vegetation (tree, shrub, herbaceous cover) (<15%)','Sparse tree (<15%)',
          'Sparse shrub (<15%)','Sparse herbaceous cover (<15%)','Tree cover, flooded, fresh or brakish water',
          'Tree cover, flooded, saline water','Shrub or herbaceous cover, flooded, fresh/saline/brakish water','Urban areas',
          'Bare areas','Consolidated bare areas','Unconsolidated bare areas','Water bodies','Permanent snow and ice']
lcols = dict(zip(lcodes,lnames))

#---------------------------------
# Read in tables
# Human influenced codes
hdc = list(range(10,41)) + [190]

# Africa
aflct = pd.read_csv(aflc)
aflct = aflct.drop(['OID_'], axis=1)
aflct['Percent'] = aflct['Count']/np.sum(aflct['Count'])*100
# Human influenced percent
afhdp = np.sum(aflct.loc[aflct['Value'].isin(hdc)].Percent)
# Urban percent
afup = aflct.loc[aflct['Value']==190].Percent.values[0]

# Asia
aslct = pd.read_csv(aslc)
aslct = aslct.drop(['OID_'], axis=1)
aslct['Percent'] = aslct['Count']/np.sum(aslct['Count'])*100
# Human influenced percent
ashdp = np.sum(aslct.loc[aslct['Value'].isin(hdc)].Percent)
# Urban percent
asup = aslct.loc[aslct['Value']==190].Percent.values[0]

# South America
salct = pd.read_csv(salc)
salct = salct.drop(['OID_'], axis=1)
salct['Percent'] = salct['Count']/np.sum(salct['Count'])*100
# Human influenced percent
sahdp = np.sum(salct.loc[salct['Value'].isin(hdc)].Percent)
# Urban percent
saup = salct.loc[salct['Value']==190].Percent.values[0]

print('Urban percent in corridors is ' + str(afup) + ' for Africa')
print('Urban percent in corridors is ' + str(asup) + ' for Asia')
print('Urban percent in corridors is ' + str(saup) + ' for South America')

sns.set(context='talk')
ax = sns.catplot(x="Value", y="Percent", data=aflct, kind='bar',height=8.27, aspect=1.25)
ax.set(xlabel='Land Cover Code', ylabel='Percent',title="Africa")

sns.set(context='talk')
ax = sns.catplot(x="Value", y="Percent", data=aslct, kind='bar',height=8.27, aspect=1.25)
ax.set(xlabel='Land Cover Code', ylabel='Percent',title="Asia")

sns.set(context='talk')
ax = sns.catplot(x="Value", y="Percent", data=salct, kind='bar',height=8.27, aspect=1.25)
ax.set(xlabel='Land Cover Code', ylabel='Percent',title="South America")

