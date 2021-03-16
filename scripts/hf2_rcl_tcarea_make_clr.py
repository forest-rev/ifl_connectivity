# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 21:26:37 2020

@author: pj276
"""
# Imports
import rasterio, np, seaborn as sns

# Only need to run this once unless changing colors
cpal = sns.color_palette("Oranges", n_colors=18)
#cpal2 = sns.color_palette("Reds", n_colors=18)

# Write clr file with one hue of orange to red for each of 18 years
rgb = []
for x, color in enumerate(cpal):
    rgb.append(['1'] + [str(np.int(np.round(i*255))) for i in color])

# Write each list element to a different file
for x, color in enumerate(rgb):
    with open('G:/My Drive/Projects/ifl_corridors/rp5k130k25k/annual_loss/tcaloss_' + str(x+1) + '.clr', 'w') as f:
        f.write('%s\n' % ' '.join(color))

#%%
# Reclassify south america rasters (32,33)
# Read raster
for n in range(1,19,1):
    rname = 'G:/My Drive/Projects/ifl_corridors/rp5k130k25k/annual_loss/sa_rtile_lybytc2000_' + str(n) + '_90m_33.tif'
    with rasterio.open(rname) as r:
        r_meta = r.profile
        a = r.read(1)
        a = np.where(a > 0, 1, 0)
    r_meta['dtype'] = 'int16'
    with rasterio.open('G:/My Drive/Projects/ifl_corridors/rp5k130k25k/annual_loss/sa_rtile_lybytc2000_' + str(n) + '_90m_33_rcl.tif', 'w', **r_meta) as dst:
        dst.write(a.astype('int16'), 1)

#%%
# Reclassify africa rasters (25,26,31,32)
# Read raster
for n in range(1,19,1):
    rname = 'G:/My Drive/Projects/ifl_corridors/rp5k130k25k/annual_loss/af_rtile_lybytc2000_' + str(n) + '_90m_33.tif'
    with rasterio.open(rname) as r:
        r_meta = r.profile
        a = r.read(1)
        a = np.where(a > 0, 1, 0)
    r_meta['dtype'] = 'int16'
    with rasterio.open('G:/My Drive/Projects/ifl_corridors/rp5k130k25k/annual_loss/af_rtile_lybytc2000_' + str(n) + '_90m_33_rcl.tif', 'w', **r_meta) as dst:
        dst.write(a.astype('int16'), 1)