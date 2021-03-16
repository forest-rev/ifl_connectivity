# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 08:58:35 2020
Calculate flux and tree cover stats every year, relative to 2000
@author: pj276
"""

# For each spatial unit (continent, country, ecoregion, thiessen polygon)
# Tabulate stats (mean, sd, and sums)

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
import geopandas as gpd, rasterio, os, ogr, sys
from shapely.geometry import box
import fiona, gdal
import pandas as pd
from rasterio.tools.mask import mask
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
# Get rastproj
rtemp = p1.rtemp #'/scratch/pj276/ifl_corridors/rp5k130k90m_neg01_as_2013_ifl/cfinal/flux/corr_2013_120_99_flux.tif'
rtemp = rasterio.open(rtemp)
rproj4 = rtemp.crs

#%%
# Read in zone shapefile
inshp = p1.inshp 
shpfile = gpd.read_file(inshp)

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

# Empty list for shapefile names
sflist = []

# Flux layer
flname = idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/flux/csum_rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl_mos.tif'

# Tree cover 2000 layer
tcname = idir + emtd + '/' + i + '_' + 'treecover2000.vrt'

# Get list of tree cover loss layers
tclosslayers = [i + '_' + 'lybytc2000_' + str(z) + '.vrt' for z in range(0,19,1)]

# Rasterized ifl name
riflname = '/scratch/pj276/ifl_corridors/ifl/ifl_2000_tropics_ssu_' + i + '_' + j + '.tif'

# Get flux raster corners
minx,miny,maxx,maxy = hfu.raster_corners(flname)
b = [box(minx, miny, maxx, maxy)]
d = {'id': [1], 'region': [i]}
dd = pd.DataFrame(data=d)
bbox = gpd.GeoDataFrame(dd, geometry=b)
bbox.crs = rproj4
bboxgeog = bbox.to_crs(shpfile.crs) 

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

# Output tiffs (otiff3 defined later)
otiff1 = idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '_flux.tif'
otiff2 = idir + emtd + '/' + j + '_' + i + '_' + k + '_' + l + '_tc2000.tif'
otiffrifl = idir + emtd + '/' + j + '_' + i + '_' + k + '_' + l + '_rifl.tif'

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
    # Check if polygon intersects with raster
    minx, miny, maxx, maxy = cx.geometry.total_bounds
    b = [box(minx, miny, maxx, maxy)]
    d = {'id': [1], 'ucode': [l]}
    dd = pd.DataFrame(data=d)
    bboxcx = gpd.GeoDataFrame(dd, geometry=b)
    bboxcx.crs = cx.crs
    ix = gpd.overlay(bboxgeog, bboxcx, how='intersection', use_sindex=True)
    if ix.shape[0] >= 1:
        # Repro shapefile
        cx = cx.to_crs(rproj4)
        # Save shapefile if needed
        if not os.path.exists(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '.shp'):
            cx.to_file(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '.shp')              
            sflist.append(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '.shp')
        # Get shapefile extent in regional sinusoidal coordinates
        minx, miny, maxx, maxy = cx.geometry.total_bounds

        # Get geometry from zone shapefile
        with fiona.open(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '.shp', "r") as shapefile:
            geoms = [feature["geometry"] for feature in shapefile]

        # Get projection info
        rwkt = gdal.Open(flname)
        rwkt = rwkt.GetProjectionRef()
        
        # Get flux extent
        fluxext = hfu.raster_corners(flname)
        # Get tree cover extent
        tcext = hfu.raster_corners(tcname)
        # Get tree cover loss extent
        tclext = hfu.raster_corners(idir + emtd + '/' + tclosslayers[0])
        # Get ifl extent
        riflext = hfu.raster_corners(riflname)
        
        # Get intersection of extents
        minx = np.max([minx, fluxext[0], tcext[0], tclext[0], riflext[0]])
        miny = np.max([miny, fluxext[1], tcext[1], tclext[1], riflext[1]])
        maxx = np.min([maxx, fluxext[2], tcext[2], tclext[2], riflext[2]])
        maxy = np.min([maxy, fluxext[3], tcext[3], tclext[3], riflext[3]])
                
        #--------------------------------------------------
        # Read in flux raster, subset by shape bounding box
        dfr = hfu.raster2array(flname, clip=[minx, miny, maxx, maxy])
        # If raster2array returns a valid output, keep processing
        if dfr != False:
            # Write clipped raster to file
            hfu.array2raster(otiff1,'#',gdal.GDT_Float32,'GTiff',dfr[0],geotrans=dfr[1],rasterproj=rwkt)
            # Use zone boundary to mask flux sum layer
            with rasterio.open(otiff1) as src:
                out_image1, out_transform = mask(src, geoms, crop=True)
                out_meta = src.meta.copy()
            # Delete rectangle cropped raster
            os.remove(otiff1)
            # Sum flux values to get total for the zone
            fluxsum = np.sum(out_image1)
            # Break out of script if there's no flux
            if fluxsum == 0:
                print('no flux, exiting script')
                sys.exit(0)
        else:
            print('invalid extent returned, exiting script')
            sys.exit(0)

        #--------------------------------------------------
        # Read in tree cover vrt, subset by shape bounding box
        dfr = hfu.raster2array(tcname, clip=[minx, miny, maxx, maxy])
        # If raster2array returns a valid output, keep processing
        if dfr != False:
            # Write clipped raster to file
            hfu.array2raster(otiff2,'#',gdal.GDT_Int16,'GTiff',dfr[0],geotrans=dfr[1],rasterproj=rwkt)
            # Use zone boundary to mask flux sum layer
            with rasterio.open(otiff2) as src:
                out_image2, out_transform = mask(src, geoms, crop=True)
                out_meta = src.meta.copy()
            # Delete rectangle cropped raster
            os.remove(otiff2)
            # Divide tree cover percent by 100 to get fractional cover,
            # then multiply by cell area (90m x 90m = 8100) to get tree covered
            # area. Divide to get km2.
            tcareasum = np.sum(out_image2/100.0*8100)/1000000
            # Subset tree cover percent to corridor area and then calculate
            # tree covered area in corridor.
            tcareaincorrsum = np.sum(out_image2[out_image1 > 0]/100.0*8100)/1000000
            # Convert to area (m2) as input to next section
            out_image2 = out_image2/100.0*8100.0
        else:
            print('invalid extent returned, exiting script')
            sys.exit(0)
        
        #--------------------------------------------------
        # Read in rasterized ifl, subset by shape bounding box
        dfr = hfu.raster2array(riflname, clip=[minx, miny, maxx, maxy])
        # If raster2array returns a valid output, keep processing
        if dfr != False:
            # Write clipped raster to file
            hfu.array2raster(otiffrifl,'#',gdal.GDT_Int16,'GTiff',dfr[0],geotrans=dfr[1],rasterproj=rwkt)
            # Use zone boundary to mask flux sum layer
            with rasterio.open(otiffrifl) as src:
                out_imagerifl, out_transform = mask(src, geoms, crop=True)
                out_meta = src.meta.copy()
            # Delete rectangle cropped raster
            os.remove(otiffrifl)
            # Subset tree cover percent to corridor area and then calculate
            # tree covered area in corridor but outside of ifls.
            tcareaincorrnoiflsum = np.sum(out_image2[(out_image1 > 0) & (out_imagerifl < 1)])/1000000.0
        else:
            print('invalid extent returned, exiting script')
            sys.exit(0)
    
#        #--------------------------------------------------
#        # Loop through annual tree cover loss layers, calculating
#        # 1. tree cover loss in the zone
#        # 2. tree cover loss in corridors in the zone
#        # 3. flux loss in the zone
#        tclosssum = []
#        tclossincorrsum = []
#        flosssum = []
#        for z, p in enumerate(tclosslayers):
#            dfr = hfu.raster2array(idir + emtd + '/' + p, clip=[minx, miny, maxx, maxy])
#            # If raster2array returns a valid output, keep processing
#            if dfr != False:
#                otiff3 = idir + emtd + '/' + j + '_' + i + '_' + k + '_' + l + '_lybtc2000_' + str(z) + '.tif'
#                # Write clipped raster to file
#                hfu.array2raster(otiff3,'#',gdal.GDT_Int16,'GTiff',dfr[0],geotrans=dfr[1],rasterproj=rwkt)
#                # Use zone boundary to mask flux sum layer
#                with rasterio.open(otiff3) as src:
#                    out_image3, out_transform = mask(src, geoms, crop=True)
#                    out_meta = src.meta.copy()
#                # Delete rectangle cropped raster
#                os.remove(otiff3)
#                # Sum tree cover loss in zone
#                tclosssum.append(np.sum(out_image3)/1000000.0)
#                # Sum tree cover loss in corridors in zone
#                tclossincorrsum.append(np.sum(out_image3[out_image1 > 0])/1000000.0)
#                # Get flux loss in zone
#                # Use out_image1 to index areas in corridors
#                # Use out_image3 to index areas of loss, which avoids dividing by 0 tree covered area in 2000
#                flosssum.append(np.sum((out_image3[out_image3 > 0.0].astype('float32')/out_image2[out_image3 > 0.0].astype('float32'))*out_image1[out_image3 > 0.0]))
#            else:
#                print('invalid extent returned, exiting script')
#                sys.exit(0)

        # One line dataframe to hold flux and other info
        #d = {'region': [i], 'scenario': [j], 'year': [k], 'ucode': [l], 'flux2000': [fluxsum], 'tca2000': [tcareasum], 'tcac2000': [tcareaincorrsum], 'tcloss': [tclosssum], 'tclossc': [tclossincorrsum], 'floss': [flosssum], 'tcacnoifl': [tcareaincorrnoiflsum]}
        d = {'region': [i], 'scenario': [j], 'year': [k], 'ucode': [l], 'tcacnoifl': [tcareaincorrnoiflsum]}
        fluxdf = pd.DataFrame(data=d)
        # Write to file
        fluxdf.to_csv(idir + '/rp5k130k90m_' + j + '_' + i + '_' + k + '_ifl/cfinal/' + fluxdir + '/' + j + '_' + i + '_' + k + '_' + l + '_lybytc2000flux_noifl.csv')

        # Delete shapefiles
        driver = ogr.GetDriverByName("ESRI Shapefile")
        for p in list(set(sflist)):
            driver.DeleteDataSource(p)