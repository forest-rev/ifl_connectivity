# -*- coding: utf-8 -*-
"""
Created on Thu Feb 08 10:58:30 2018
This script calculates regularly spaced points for each polygon in a shapefile.
It is the first of three parts. It maps the points.
The script does check for existing point shapefiles, skipping them if they exist.
@author: pj276
"""

#%%
import sys, os, fiona, ogr, glob
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
sys.path.append('/home/pj276/projects/ifl_corridors/code/functions/')
import hf_utilities as hfu
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
#import p1_10k_100k_041118 as p1
#import p1_20k_100k_052518 as p1
#import p1_10k_100k_060418 as p1
#import p1_10k_300k_060518 as p1
#import p1_90m_110718_as p1
#import p1_90m_ifl_2000_sa as p1
import param_file as p1

#%%
# Set up inputs
odir = p1.odir
pointdir = p1.pointdir
#'/projects/above_gedi/pjantz/ifl_corridors_data/outputs/regular_points/rp10ks1'
xspac = p1.xspac
#10000 # x spacing in meters
yspac = p1.yspac
#10000 # y spacing in meters
#xsnap = -3662648.020845168270171 # south america custom sinusoidal upper left x
#ysnap = 3320113.397940382361412 # south america custom sinusoidal upper left y
fieldname = p1.fieldname
#'GRIDCODE'
tdist = p1.tdist
#300000 # maximum point distance to consider, in meters
#inpoly = '/projects/above_gedi/pjantz/ifl_corridors_data/hifo/hland_tropics_ssu_sa.shp'
inpoly = p1.inpoly
#'/projects/above_gedi/pjantz/ifl_corridors_data/hifo/hland_tropics_ssu_sa_subset.shp'
dname = p1.dname
#'10k' # string corresponding to x,y spacing for naming
csvout = p1.csvout
#'/projects/above_gedi/pjantz/ifl_corridors_data/outputs/regular_points/rp10ks1/ptnbs_10k.csv'
pointlog = p1.pointlog

#%%
# Get projection
t1 = fiona.open(inpoly)
myproj = t1.crs
t1.close()
t1 = None

#%%
# Get unique gridcodes and filter out 0 gridcodes
# Empty list
arealist = []
reclist = []
# Read in polygon
with fiona.open(inpoly) as input:
    for rec in input:
        reclist.append(rec['properties'][fieldname])
        arealist.append(rec['properties']['POLY_AREA'])
idrec = [reclist,arealist]
df = pd.DataFrame(idrec)
df = df.transpose()
df.columns = [fieldname, 'AREA']
df = df[df[fieldname] > 0]
dfs = df.groupby(fieldname).sum()
# Get unique gridcode
ugc = list(dfs.index.unique())

#%%
##pcxlist = []
##pcylist = []
##gclist = []
##ptidlist = []
# Loop throught polygons, create regular coordinates, and append to list
for g in ugc:
    # Convert to integer
    g = int(g)
    print(g)
    # Make directory to hold points
    if not os.path.exists(odir + pointdir + os.sep + str(g)):
        os.mkdir(odir + pointdir + os.sep + str(g))
    # Intermediate shapefile name
    dissoshape = odir + pointdir + os.sep + str(g) + os.sep + os.path.basename(inpoly).split('.')[0] + '_' + str(g) + '.shp'
    # Final point shapefile name
    opath = dissoshape.replace('.shp','_' + dname + '.shp')
    # If all 5 shapfile files aren't there, delete what is
    # If they are all there, the next step should skip
    spat = opath.replace('.shp','*')
    nfiles = glob.glob(spat)
    if len(nfiles) > 0 and len(nfiles) < 5:
        for f in nfiles:
            os.remove(f)
    # If the final point file exists, skip it
    if not os.path.exists(opath):
        # If final point shapefile doesn't exist already, create it.
        # Use fiona and geopandas to dissolve features if there are multiple ones with the same gridcode
        reader = fiona.open(inpoly)
        xx = gpd.GeoDataFrame.from_features((x for x in reader if x['properties'][fieldname]==g))
        xx = xx.dissolve(by=fieldname, aggfunc='sum')
        xx.crs = myproj
        xx[fieldname] = g
        # Write polygon to shapefile
        xx.to_file(dissoshape)
        reader.close()
        # Create regular points array
        #pco = hfu.regular_points(dissoshape, fieldname, g, xspac, yspac, dosnap=1, snap_x=xsnap, snap_y=ysnap)
        pco = hfu.regular_points_from_center(dissoshape, fieldname, g, xspac, yspac)
        # Convert pco to gpd object
        df = pd.DataFrame({'uid': range(0,len(pco)), 'xco': [z[0] for z in pco], 'yco': [z[1] for z in pco]})
        df['Coordinates'] = pco
        df['Coordinates'] = df['Coordinates'].apply(Point)
        gdf_point = gpd.GeoDataFrame(df, geometry='Coordinates')
        gdf_poly = gpd.read_file(dissoshape)
        gdf_point.crs = gdf_poly.crs
        ajoin = gpd.tools.sjoin(gdf_point, gdf_poly, how='left', op='within')
        ajoin = ajoin.groupby('index_right')
        ajoin = list(ajoin)[0][1]
        ajoin = pd.DataFrame(ajoin)
        ajoin['id'] = range(0,len(ajoin.uid))
        ajoin = ajoin[['id','xco', 'yco',fieldname]]
        ajoin['Coordinates'] = zip(ajoin.xco.tolist(),ajoin.yco.tolist())
        ajoin['Coordinates'] = ajoin['Coordinates'].apply(Point)
        ajoin = gpd.GeoDataFrame(ajoin, geometry='Coordinates')
        ajoin.crs = gdf_poly.crs
        ajoin.rename(columns={fieldname:'GCODE'}, inplace=True)
        if len(ajoin.id) < 3:
            xspac = 4000
            yspac = 4000
            pco = hfu.regular_points_from_center(dissoshape, fieldname, g, xspac, yspac)
            # Convert pco to gpd object
            df = pd.DataFrame({'uid': range(0,len(pco)), 'xco': [z[0] for z in pco], 'yco': [z[1] for z in pco]})
            df['Coordinates'] = pco
            df['Coordinates'] = df['Coordinates'].apply(Point)
            gdf_point = gpd.GeoDataFrame(df, geometry='Coordinates')
            gdf_poly = gpd.read_file(dissoshape)
            gdf_point.crs = gdf_poly.crs
            ajoin = gpd.tools.sjoin(gdf_point, gdf_poly, how='left', op='within')
            ajoin = ajoin.groupby('index_right')
            ajoin = list(ajoin)[0][1]
            ajoin = pd.DataFrame(ajoin)
            ajoin['id'] = range(0,len(ajoin.uid))
            ajoin = ajoin[['id','xco', 'yco',fieldname]]
            ajoin['Coordinates'] = zip(ajoin.xco.tolist(),ajoin.yco.tolist())
            ajoin['Coordinates'] = ajoin['Coordinates'].apply(Point)
            ajoin = gpd.GeoDataFrame(ajoin, geometry='Coordinates')
            ajoin.crs = gdf_poly.crs
            ajoin.rename(columns={fieldname:'GCODE'}, inplace=True)
            # Write file paths to text file
            with open(pointlog, 'a') as afile:
                afile.write(dissoshape + 'at 4k \n')
            if len(ajoin.id) < 3:
                xspac = 3000
                yspac = 3000
                pco = hfu.regular_points_from_center(dissoshape, fieldname, g, xspac, yspac)
                # Convert pco to gpd object
                df = pd.DataFrame({'uid': range(0,len(pco)), 'xco': [z[0] for z in pco], 'yco': [z[1] for z in pco]})
                df['Coordinates'] = pco
                df['Coordinates'] = df['Coordinates'].apply(Point)
                gdf_point = gpd.GeoDataFrame(df, geometry='Coordinates')
                gdf_poly = gpd.read_file(dissoshape)
                gdf_point.crs = gdf_poly.crs
                ajoin = gpd.tools.sjoin(gdf_point, gdf_poly, how='left', op='within')
                ajoin = ajoin.groupby('index_right')
                ajoin = list(ajoin)[0][1]
                ajoin = pd.DataFrame(ajoin)
                ajoin['id'] = range(0,len(ajoin.uid))
                ajoin = ajoin[['id','xco', 'yco',fieldname]]
                ajoin['Coordinates'] = zip(ajoin.xco.tolist(),ajoin.yco.tolist())
                ajoin['Coordinates'] = ajoin['Coordinates'].apply(Point)
                ajoin = gpd.GeoDataFrame(ajoin, geometry='Coordinates')
                ajoin.crs = gdf_poly.crs
                ajoin.rename(columns={fieldname:'GCODE'}, inplace=True)
                with open(pointlog, 'a') as afile:
                    afile.write(dissoshape + 'at 3k \n')
        # Get projection info from input shapefile
##        f = open(os.path.dirname(dissoshape) + os.sep + os.path.basename(dissoshape).replace('shp','prj'))
##        fr = f.read() # Read wkt to variable
##        f = None # Close file
##        srs=osr.SpatialReference(wkt=fr) # convert to object
        ajoin.to_file(opath)
        # Write coordinates to shapefile (the existence check isn't necessary)
##        if os.path.exists(opath):
##            shpDriver = ogr.GetDriverByName('ESRI Shapefile')
##            shpDriver.DeleteDataSource(opath)
##            shpDriver = None
##        hfu.coordinates2point(pco,opath,srs)
        # Add gridcode as field
##        hfu.add_shp_field(opath, 'GCODE', 'N', 8, int(g), opath)
        # Delete merged polygon
        if os.path.exists(dissoshape):
            shpDriver = ogr.GetDriverByName('ESRI Shapefile')
            shpDriver.DeleteDataSource(dissoshape)
            shpDriver = None
        # Append gridcode and coordinates to list
        ##pcxlist.append(ajoin.xco.tolist())
        ##pcylist.append(ajoin.yco.tolist())
        ##gclist.append([g]*len(ajoin.id))
        ##ptidlist.append(range(0,len(ajoin.id)))
        #pco = zip(pco, [g]*len(pco))
        #pclist.append(pco)
        adict1 = {'xco': ajoin.xco.tolist(), 'yco': ajoin.yco.tolist(), 'gcode': [g]*len(ajoin.id), 'ptid': range(0,len(ajoin.id))}
        adf1 = pd.DataFrame(adict1)
        adf1.index.name = 'rid'
        csvout1 = os.path.dirname(csvout) + '/' + os.path.basename(csvout).split('.')[0] + '_' + str(g) + '_temp.csv'
        adf1.to_csv(csvout1)
        
# If any files were processed, put the info into a csv
##if len(pcxlist) > 0:
##    pcxlist = list(chain.from_iterable(pcxlist))
##    pcylist = list(chain.from_iterable(pcylist))
##    gclist = list(chain.from_iterable(gclist))
##    ptidlist = list(chain.from_iterable(ptidlist))
##    adict = {'xco': pcxlist, 'yco': pcylist, 'gcode': gclist, 'ptid': ptidlist}
##    adf = pd.DataFrame(adict)
##    adf.index.name = 'rid'
##    adf.to_csv(csvout)

