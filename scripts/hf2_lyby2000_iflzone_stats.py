# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 16:25:44 2020
Calculate flux and tree cover stats every year, relative to 2000
@author: pj276
"""

# For each spatial unit (continent, country, ecoregion, thiessen polygon)
# Tabulate stats

# 1.a. read or make lossyear vrt from 1.6 gfc data
    # b. read or make year 2000 flux vrt
# 2.a. tabulate year 2000 flux
    # b. tabulate year 2000 forest cover
# 3.loop through years
    # a. extract loss patches
    # b. overlay with year 2000 forest cover and tabulate 
    # c. overlay with year 2000 forest cover in corridors and tabulate
    # d. overlay with year 2000 flux and tabulate

#%%
# Imports
import sys
import geopandas as gpd
import pandas as pd
import numpy as np
sys.path.append('/home/pj276/projects/ifl_corridors/code/functions/')
import hf_utilities as hfu
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
import param_file as p1

#%%
# Top level dir
idir = p1.idir #'/scratch/pj276/ifl_corridors'

# Directory that will hold output direct flux files
fluxdir = p1.fluxdir

# Edge matched tile directory for gfc data
emtd = p1.emtd

#%%
# Read in zone shapefile
inshp = p1.inshp 
shpfile = gpd.read_file(inshp,  encoding="utf-8")

# Specify field
field = p1.field

# Get list of unique names/codes
ucode = np.unique(shpfile[field].tolist())

# Job array index
jai = int(sys.argv[1]) - 1 # subtract for python zero indexing
    
#%%
# Region
i = p1.region #['af','as','sa']
# Scenario
j = p1.scenario #['neg2', 'neg01', 'pos2']

# Year
k = '2000'

# Zone raster
zr = '/scratch/pj276/ifl_corridors/iflzones/' + i + '_rtile_cs' + k + j + '_90m_calloc.tif'

# Euclidean zone raster
ezr = '/scratch/pj276/ifl_corridors/iflzones/' + i + '_rtile_' + k  + '_90m_ealloc.tif'

# Flux layer
flname = idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/flux/csum_rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl_mos.tif'

# Tree cover 2000 layer
tcname = idir + emtd + '/' + i + '_' + 'treecover2000.vrt'

# Datamask 2000 layer
dmname = '/scratch/pj276/ifl_corridors/gfc/gfc1_0/gfc1_0_' + i + '_ss_emt_90m/mos_datamask_' + k + '.vrt'

# Get list of tree cover loss layers
tclosslayers = [i + '_' + 'lybytc2000_' + str(z) + '.vrt' for z in range(0,19,1)]

# Test if unique code ucode is string or unicode, treat differently if so
if isinstance(ucode[jai], str):
    l = ucode[jai]
    lorig = l
    if ' ' in l:
        l = l.replace(' ', '_')
elif isinstance(ucode[jai], unicode):
    l = ucode[jai]
    lorig = l
    if ' ' in l:
        l = l.replace(' ', '_')
else:
    l = ucode[jai]
    lorig = l
    l = np.str(np.int(l))

#%%
# Get shapes
cx = shpfile.loc[[s for s, n in enumerate(shpfile[field]) if n == lorig],]
# If at least one shape is returned, keep processing
if cx.shape[0] > 0:
    # Get target column and geometry column
    cx = cx[[field,'geometry']]
    # Dissolve if there are multiple shapes
    if cx.shape[0] > 1:
        cx = cx.dissolve(by=field, aggfunc='first')
        cx[field] = cx.index
    # Fix simple self intersections if necessary
    if cx.is_valid.bool() == False:
        cx = cx.buffer(0)
    # Get geometry bounds and buffer by one cell
    minx, miny, maxx, maxy = cx.geometry.total_bounds
    minx = minx-90.0
    miny = miny-90.0
    maxx = maxx+90.0
    maxy = maxy+90.0
           
    # Get flux extent
    fluxext = hfu.raster_corners(flname)
    # Get tree cover extent
    tcext = hfu.raster_corners(tcname)
    # Get tree cover loss extent
    tclext = hfu.raster_corners(idir + emtd + '/' + tclosslayers[0])
    
    # Get intersection of extents
    minx = np.max([minx, fluxext[0], tcext[0], tclext[0]])
    miny = np.max([miny, fluxext[1], tcext[1], tclext[1]])
    maxx = np.min([maxx, fluxext[2], tcext[2], tclext[2]])
    maxy = np.min([maxy, fluxext[3], tcext[3], tclext[3]])

    # Read in zone raster subset
    zrast = hfu.raster2array(zr, clip=[minx, miny, maxx, maxy])
    
    # Read in euclidean zone raster subset
    ezrast = hfu.raster2array(ezr, clip=[minx, miny, maxx, maxy])
    
    # Read in flux raster, subset by zone bounding box
    frast = hfu.raster2array(flname, clip=[minx, miny, maxx, maxy])
    # If raster2array returns a valid output, keep processing
    if frast != False:
        zrast = zrast[0]
        frast = frast[0]
        # Sum flux values to get total for the zone
        fluxsum = np.sum(frast[(zrast==lorig) & (ezrast>0)])
        # Break out of script if there's no flux
        if fluxsum == 0:
            print('no flux, exiting script')
            sys.exit(0)
    else:
        print('invalid extent returned, exiting script')
        sys.exit(0)

    #--------------------------------------------------
    # Read in tree cover vrt, subset by shape bounding box
    tcrast = hfu.raster2array(tcname, clip=[minx, miny, maxx, maxy])
    # If raster2array returns a valid output, keep processing
    if tcrast != False:
        tcrast = tcrast[0]
        # Divide tree cover percent by 100 to get fractional cover,
        # then multiply by cell area (90m x 90m = 8100) to get tree covered
        # area. Divide to get km2.
        tcareasum = np.sum(tcrast[(zrast==lorig) & (ezrast>0)]/100.0*8100)/1000000
        # Subset tree cover percent to corridor area and then calculate
        # tree covered area in corridor.
        tcareaincorrsum = np.sum(tcrast[(frast > 0) & (zrast==lorig) & (ezrast>0)]/100.0*8100)/1000000
        # Convert to area (m2) as input to next section
        tcrast = tcrast/100.0*8100.0
    else:
        print('invalid extent returned, exiting script')
        sys.exit(0)
    
    #--------------------------------------------------
    # Read in datamask vrt, subset by shape bounding box
    dmrast = hfu.raster2array(dmname, clip=[minx, miny, maxx, maxy])
    # If raster2array returns a valid output, keep processing
    if dmrast != False:
        dmrast = dmrast[0]
        # Get land area in zone
        dmareasum = np.sum(dmrast[(zrast==lorig) & (ezrast>0) & (dmrast==1)]*8100.0)/1000000.0
        # Subset tree cover percent to corridor area and then calculate
        # tree covered area in corridor.
        dmareaincorrsum = np.sum(dmrast[(frast > 0) & (zrast==lorig) & (ezrast>0) & (dmrast==1)]*8100)/1000000.0
        del dmrast
    else:
        print('invalid extent returned, exiting script')
        sys.exit(0)
    
    #--------------------------------------------------
    # Loop through annual tree cover loss layers, calculating
    # 1. tree cover loss in the zone
    # 2. tree cover loss in corridors in the zone
    # 3. flux loss in the zone
    tclosssum = []
    tclossincorrsum = []
    flosssum = []
    for z, p in enumerate(tclosslayers):
        tclrast = hfu.raster2array(idir + emtd + '/' + p, clip=[minx, miny, maxx, maxy])
        # If raster2array returns a valid output, keep processing
        if tclrast != False:
            tclrast = tclrast[0]
            # Sum tree cover loss in zone
            tclosssum.append(np.sum(tclrast[(zrast==lorig) & (ezrast>0)])/1000000.0)
            # Sum tree cover loss in corridors in zone
            tclossincorrsum.append(np.sum(tclrast[(frast > 0) & (zrast==lorig) & (ezrast>0)])/1000000.0)
            # Get flux loss in zone
            # Use out_image1 to index areas in corridors
            # Use out_image3 to index areas of loss, which avoids dividing by 0 tree covered area in 2000
            flosssum.append(np.sum((tclrast[(tclrast > 0.0) & (zrast==lorig) & (ezrast>0)].astype('float32')/tcrast[(tclrast > 0.0) & (zrast==lorig)].astype('float32'))*frast[(tclrast > 0.0) & (zrast==lorig) & (ezrast>0)]))
        else:
            print('invalid extent returned, exiting script')
            sys.exit(0)
    # One line dataframe to hold flux and other info
    d = {'region': [i], 'scenario': [j], 'year': [k], 'ucode': [l], 'flux2000': [fluxsum], 'tca2000': [tcareasum], 'tcac2000': [tcareaincorrsum], 'tcloss': [tclosssum], 'tclossc': [tclossincorrsum], 'floss': [flosssum], 'landarea': [dmareasum], 'landareac': [dmareaincorrsum]}
    fluxdf = pd.DataFrame(data=d)
    # Write to file
    fluxdf.to_csv(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '_lybytc2000flux.csv')
    
