# -*- coding: utf-8 -*-
"""
Created on Sun Jul 02 20:40:58 2017
This script takes a fishnet shapefile, an index shapefile of existing rasters
created by gdaltindex, and mosaics the rasters together using a maximum
operator within the fishnet polygons. The resulting max-mosaiced rasters
should edge match.
@author: pj276
"""
import shapefile as shp, geopandas as gpd
import sys, os, subprocess, gdal, rasterio
from shapely.geometry import box
sys.path.append('/home/pj276/projects/ifl_corridors/code/functions/')
# For hf_utilities, the modules are available within the functions
# when they're used. Otherwise, one needs to prefix with hfu.
# e.g. hfu.shp.Writer
import hf_utilities as hfu
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
sys.path.append('/home/pj276/projects/undp_connectivity/code/parameter_files/')
#import p_102418 as p1
import param_file as p1

###############################################################################
# Set up arguments
foresttiles = p1.foresttiles # List of strings corresponding to forest layers of interest
continent = p1.continent # Target continent
fnet = p1.fnet # Fishnet

###############################################################################

# Loop through layers and edge match
for k in foresttiles:
    # Set up variables for layer
    fmosindex = p1.tld + p1.pfdd + '/tind_' + k + '_' + continent + '.shp' # Tile index describing original raster footprints
    txtlistpre = p1.tld + p1.pfdd + '/txtlist_' + k + 'A_' + continent + '_' # this is a text file
    vrt2pre = p1.tld + p1.pfdd + '/vrtshingle' + k + '_' + continent + '_' # this is a vrt
    rtilepre = p1.tld + p1.emtd + '/' + continent + '_rtile_' + k.strip('_') + '_' # this is a tiff
    mosvrtlog = p1.tld + p1.emtd + '/hf_mosaic_vrt_log_' + k + '_' + continent + '.txt' # this is a text file    

    # Make tile index if it doesn't already exist
    if not os.path.exists(fmosindex):
        # Translate tindex file name to raster file name search pattern
        df = gpd.GeoDataFrame(columns=['location','geometry'])
        # List files and get those that match search pattern
        flist = os.listdir(os.path.dirname(fmosindex))
        flist = [os.path.dirname(fmosindex) + '/' + i for i in flist if k in i and i.endswith('.tif')]
        if len(flist) == 0:
              sys.exit('no files to process')
        for fname in flist:
            bounds = rasterio.open(fname).bounds
            df = df.append({'location':fname, 'geometry': box(bounds[0], bounds[1], bounds[2], bounds[3])},ignore_index=True)
        # Use last raster file to get crs and give it to shapefile
        df.crs = rasterio.open(fname).crs
        # Check for file again before writing
        if not os.path.exists(fmosindex):
            df.to_file(fmosindex)
    
    # Edge match tiles
    for j in range(0,hfu.getfc(fnet)):
        rtile = rtilepre + str(j) + ".tif"
        if not os.path.exists(rtile):         
            # List to hold names of rasters that overlap fishnet tile
            ilist = []
            for i in range(0,hfu.getfc(fmosindex)):
                ilist.append(hfu.intersect(fnet,fmosindex,j,i))
            ilist = [x for x in ilist if x is not None]
            if len(ilist) > 0:        
                # Write to file for input to gdalbuildvrt
                #txtlist = "/scratch/pj276/ifl_corridors/geodata/hires_global_v1_0_2013/txtlistA_" + str(j) + ".txt"
                #txtlist = "/scratch/pj276/ifl_corridors/geodata/WaterMask2010_UMD/txtlistA_" + str(j) + ".txt"
                txtlist = txtlistpre + str(j) + ".txt"
                with open(txtlist,"a") as tf:
                    for li in ilist:
                        tf.write(li + "\n")    
                # Get polygon coordinate mins and maxes
                sf = shp.Reader(fnet)
                ashape = sf.shapes()[j]
                xpts = [] # list to hold x coords
                ypts = [] # list to hold y coords
                # Populate lists of x and y coords
                for u in ashape.points:
                    xpts.append(u[0])
                    ypts.append(u[1])
                # Get mins and maxes
                xmin = str(min(xpts))
                xmax = str(max(xpts))
                ymin = str(min(ypts))
                ymax = str(max(ypts))
                
                # Build a shingled vrt and then subset with read as array for 
                # to extract edge matched rasters based on fishnet tiles. This
                # might work for the entire domain but test with two rasters to start with.
                #vrt2 = "/scratch/pj276/ifl_corridors/geodata/WaterMask2010_UMD/vrtshingle" + str(j) + ".vrt"
                vrt2 = vrt2pre + str(j) + ".vrt"
                subprocess.call(["gdalbuildvrt", "-separate", "-te", xmin, ymin, xmax, ymax, "-input_file_list", txtlist, vrt2])
                
                # Take the max of rasters where they overlap within the extent of the fishnet tile
                #rtile = "/scratch/pj276/ifl_corridors/geodata/WaterMask2010_UMD/sa_rtile_" + str(j) + ".tif"
                hfu.tile_shingle(vrt2, rtile, p1.gdpt, p1.nppt, 'max')
                # Delete vrt files
            #    os.remove(vrt1)
                os.remove(vrt2)
                os.remove(txtlist)
                #with open("/scratch/pj276/ifl_corridors/geodata/WaterMask2010_UMD/hf_mosaic_vrt_log.txt", "a") as text_file:
                with open(mosvrtlog, "a") as text_file:
                    text_file.write("Finished tile " + str(j) + " of " + str(hfu.getfc(fnet)) + "\n")
            else:
                #with open("/scratch/pj276/ifl_corridors/geodata/WaterMask2010_UMD/hf_mosaic_vrt_log.txt", "a") as text_file:
                with open(mosvrtlog, "a") as text_file:
                    text_file.write("Tile " + str(j) + " of " + str(hfu.getfc(fnet)) + " is empty" + "\n")

