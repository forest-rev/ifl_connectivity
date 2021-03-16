# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 20:18:16 2018
This script coarsens rasters
@author: pj276
"""
import sys, os, subprocess, gdal, numpy as np
sys.path.append('/home/pj276/projects/ifl_corridors/code/functions/')
import hf_utilities as hfu
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
#import p1_90m_110718_as as p1
#import p1_90m_ifl_2000_sa as p1
#import p1_90m_ifl_2013_sa as p1
#import param_file as p1
#import flux540neg01af as p1
import param_file as p1

# Input dir
#idir = '/scratch/pj276/gfc/gfc1_4_as_ss_emt'
#idir = '/scratch/pj276/gfc/gfc1_4_as_ss_cs'
#idir = '/scratch/pj276/rp5k100k90m_' + ct + '_' + ayear + '_ifl/cfinal/flux'
#idir = '/scratch/pj276/gfc/gfc1_4_' + ct + '_ss_cs_neg2'
aggdir = p1.aggdir

# Output dir
#odir = '/scratch/pj276/gfc/gfc1_4_as_ss_emt_90m'
#odir = '/scratch/pj276/gfc/gfc1_4_as_ss_cs_90m'
#odir = '/scratch/pj276/rp5k100k90m_' + ct + '_' + ayear + '_ifl/cfinal/flux'
#odir = '/scratch/pj276/gfc/gfc1_4_' + ct + '_ss_cs_neg2_90m'
aggodir = p1.aggodir

if not os.path.exists(aggodir):
    os.mkdir(aggodir)

# The amount by which to expand bounding box dimensions
#exvec = [-90, -90, 90, 90]
exvec = p1.exvec #[-540, -540, 540, 540]
# Cell dim
#cdim = "90"
aggdim = p1.aggdim #"540"

# No data val
ndval = p1.ndval #"255"
    
# Source file name pattern
#oldname = 'treecover2013'
#oldname = 'cs2000neg01'
#oldname = 'cs2013neg01'
#srcaggname = 'rp5k100k90m_' + ct + '_' + ayear + '_ifl_mos'
#oldname = 'cs' + ayear + 'neg2'
srcaggname = p1.srcaggname

# Replacement file name pattern
#rname = 'treecover2013_90m'
#rname = 'cs2000neg01_90m'
#rname = 'cs2013neg01_90m'
#dstaggname = 'rp5k100k540m_' + ct + '_' + ayear + '_ifl_mos'
#rname = 'cs' + ayear + 'neg2_90m'
dstaggname = p1.dstaggname

# Loop through loss years
for z, p in enumerate(srcaggname):
    # List files (making sure to return full paths)
    flist = os.listdir(aggdir)
    flist = [aggdir + '/' + i for i in flist if p in i and i.endswith('.tif')]
    # Loop through files
    for k, ifile in enumerate(flist):
        oraster = aggodir + '/' + os.path.basename(ifile).replace(p, dstaggname[z])
        if not os.path.exists(oraster):
            # xmin, ymin, xmax, ymax
            bbox = hfu.raster_corners(ifile) # get input file bounding box
            newbb = [] # set up new bounding box dimensions
            for i,j in enumerate(bbox):
                newbb.append(bbox[i] + exvec[i])
            
            # Build vrt with buffer
            ovrt = aggodir + '/tempvrt_' + p + '_' + str(k) + '.vrt'
            
            # Build a vrt from a file list
            acmd = ["gdalbuildvrt", "-te", str(newbb[0]), str(newbb[1]), str(newbb[2]), str(newbb[3]), ovrt]
            for myfile in flist:
                acmd.append(myfile)
            subprocess.call(acmd)
            
            # If aggregating datamask files, reclassify water values (2) to 0
            if 'datamask' in p:
                # Read vrt to array
                dm = hfu.raster2array(ovrt)[0]
                # Reclassify
                dm = np.where(dm == 2, 0, dm)
                # Name for temporary tiff
                ttiff = aggodir + '/ttiff_' + p + '_' + str(k) + '.tif'
                # Write to tiff
                hfu.array2raster(ttiff, ovrt, gdal.GDT_Byte, 'GTiff', dm)
                # Aggregate to new resolution using average operator
                ofile = aggodir + '/temp_90m_' + p + '_' + str(k) + '.tif'
                subprocess.call(["gdalwarp", "-multi", "-wo", "NUM_THREADS=ALL_CPUS", "-tr", aggdim, aggdim, "-r", "average", "-srcnodata", ndval, "-ot", "Float32", ttiff, ofile])
                # Create vrt from aggregated file, clipping off the collar
                aggvrt = aggodir + '/agg_temp_' + p + str(k) + '.vrt'
                subprocess.call(["gdalbuildvrt", "-te", str(bbox[0]), str(bbox[1]), str(bbox[2]), str(bbox[3]), aggvrt, ofile])
                # Read raster 2 array
                aggfile = hfu.raster2array(aggvrt)[0]
                # Read array 2 raster
                hfu.array2raster(oraster, aggvrt, gdal.GDT_Float32, "GTiff", aggfile)
                # Clean up 
                os.remove(aggvrt)
                os.remove(ofile)
                os.remove(ovrt)
                os.remove(ttiff)
            elif 'lybytc2000' in p:
                # If aggregating lossyear by treecover 2000 layers, take the average
                # and then multiply by 9 to get the sum
                ofile = aggodir + '/temp_90m_' + p + '_' + str(k) + '.tif'
                subprocess.call(["gdalwarp", "-multi", "-wo", "NUM_THREADS=ALL_CPUS", "-tr", aggdim, aggdim, "-r", "average", "-srcnodata", ndval, "-ot", "Float32", ovrt, ofile])
                # Create vrt from aggregated file, clipping off the collar
                aggvrt = aggodir + '/agg_temp_' + p + str(k) + '.vrt'
                subprocess.call(["gdalbuildvrt", "-te", str(bbox[0]), str(bbox[1]), str(bbox[2]), str(bbox[3]), aggvrt, ofile])
                # Read raster 2 array
                aggfile = hfu.raster2array(aggvrt)[0]
                aggfile = aggfile*9
                aggfile = (aggfile+0.5).astype('int16')
                # Read array 2 raster
                hfu.array2raster(oraster, aggvrt, gdal.GDT_Int16, "GTiff", aggfile)
                # Clean up 
                os.remove(aggvrt)
                os.remove(ofile)
                os.remove(ovrt)
            else:
                # Aggregate to new resolution using average operator
                ofile = aggodir + '/temp_90m_' + p + '_' + str(k) + '.tif'
                subprocess.call(["gdalwarp", "-multi", "-wo", "NUM_THREADS=ALL_CPUS", "-tr", aggdim, aggdim, "-r", "average", "-srcnodata", ndval, "-ot", "Float32", ovrt, ofile])
                # Create vrt from aggregated file, clipping off the collar
                aggvrt = aggodir + '/agg_temp_' + p + str(k) + '.vrt'
                subprocess.call(["gdalbuildvrt", "-te", str(bbox[0]), str(bbox[1]), str(bbox[2]), str(bbox[3]), aggvrt, ofile])
                # Read raster 2 array
                aggfile = hfu.raster2array(aggvrt)[0]
                # Read array 2 raster
                hfu.array2raster(oraster, aggvrt, gdal.GDT_Float32, "GTiff", aggfile)
                # Clean up 
                os.remove(aggvrt)
                os.remove(ofile)
                os.remove(ovrt)
