# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 19:58:51 2017
Reproject forest cover.
Can run gdaltindex -write_absolute_path for2013_index.shp *ss.tif after creating
layers to get footprints. Script assumes lat/long input projection.
@author: pj276
"""
import os, sys, math, subprocess
import gdal, gdalconst
sys.path.append('/home/pj276/projects/ifl_corridors/code/functions/')
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
sys.path.append('/home/pj276/projects/undp_connectivity/code/parameter_files/')
import param_file as p1

###############################################################################
# Set up arguments 
# Text file holding file names
#foresttiles = '/scratch/pj276/ifl_corridors/geodata/hires_global_v1_0/foresttiles.txt'
foresttiles = p1.foresttiles
# Output directory for projected forest data
#odfor = '/scratch/pj276/ifl_corridors/geodata/hires_global_v1_0'
odfor = p1.odfor
# File for reference projection
#prjref = '/scratch/pj276/ifl_corridors/geodata/hinterland_2013/hland_tropics_ss.prj'
#prjref = p1.prjref
prjrast = p1.prjrast

# Python log file
#pylogfor = "/scratch/pj276/ifl_corridors/geodata/hires_global_v1_0/pylog.txt"
#pylogfor = p1.pylogfor

# Set snap coordinates
#snap_y = 3320131.6297796257 + (10000*30)# expand ymax snap coordinate a bit
snap_y = p1.snap_y
#snap_x = -10464713.986991046 - (10000*30)# expand xmin snap coordinate a bit
snap_x = p1.snap_x
# Cell dimensions
xint = p1.xint
yint = p1.yint
# Projection suffix
#projsuffix = "_ss.tif"
projsuffix = p1.projsuffix

# Continent
continent = p1.continent

###############################################################################
# Loop through treecover, datamask, gain, and lossyear files
for k in foresttiles:
    
    # List year 2000 tree cover files (making sure to return full paths)
    fl = os.listdir(p1.tld + p1.fcd)
    fl = [p1.tld + p1.fcd + '/' + i for i in fl if k in i and i.endswith('.tif')]
    
    # Loop through file tiles and reproject
    for j in fl:
        # Reproject and resample layers. 
        src_filename = j
        src = gdal.Open(src_filename, gdalconst.GA_ReadOnly)
        src_proj = src.GetProjection()
        src_geotrans = src.GetGeoTransform()
        
        # Create output file name
        dst_filename = odfor + os.sep + os.path.basename(src_filename).replace('.tif',projsuffix)
        
        # Run if output doesn't exist
        if not os.path.exists(dst_filename):    
            # If xmax of original raster is > 180, clip it to avoid wraparound
            # Same with xmin if it is < -180.
            # Need to get xmax from gdalinfo since I get wraparound when I try to calculate from geotransform
            gi = gdal.Info(src_filename).split("\n")
            gimax = [i for i in gi if "Upper Right" in i]
            #gi = float(gi[0].split("(")[1].split(" ")[1].split(",")[0])
            gimax = float(gimax[0].split("(")[1].replace(" ","").split(",")[0])
            # Get xmin
            gimin = [i for i in gi if "Upper Left" in i]
            #gi = float(gi[0].split("(")[1].split(" ")[1].split(",")[0])
            gimin = float(gimin[0].split("(")[1].replace(" ","").split(",")[0])
            # Clip rasters if they're beyond domain edges
            if gimax > 180:
                # Specify new xmax 
                newcoord = "179.99"
                newof = os.path.dirname(src_filename) + os.sep + os.path.basename(src_filename).split(".")[0] + "_v2.tif"
                # Subset raster. Note subset window given in ulx, uly, lrx, lry order
                subprocess.call(["gdal_translate", "-co", "COMPRESS=LZW", "-projwin", str(src_geotrans[0]), str(src_geotrans[3]), newcoord, str(src_geotrans[3] + src_geotrans[5]*src.RasterYSize), src_filename, newof])
                # Reproject and resample layers. 
                src = gdal.Open(newof, gdalconst.GA_ReadOnly)
                src_proj = src.GetProjection()
                src_geotrans = src.GetGeoTransform()
            if gimin < -180:
                # Specify new xmin
                newcoord = "-179.99"
                newof = os.path.dirname(src_filename) + os.sep + os.path.basename(src_filename).split(".")[0] + "_v2.tif"
                # Subset raster. Note subset window given in ulx, uly, lrx, lry order
                subprocess.call(["gdal_translate", "-co", "COMPRESS=LZW", "-projwin", newcoord, str(src_geotrans[3]), str(src_geotrans[0] + src_geotrans[1]*src.RasterXSize), str(src_geotrans[3] + src_geotrans[5]*src.RasterYSize), src_filename, newof])
                # Reproject and resample layers. 
                src = gdal.Open(newof, gdalconst.GA_ReadOnly)
                src_proj = src.GetProjection()
                src_geotrans = src.GetGeoTransform()
        
            # Get default extent of virtually projected raster
            # Get projection from template raster
            rprjobj = gdal.Open(prjrast, gdalconst.GA_ReadOnly)
            rprj_wkt = rprjobj.GetProjection() 
            
            # Virtual VRT
            tfile = gdal.AutoCreateWarpedVRT(src, src_proj, rprj_wkt)
            dst_xsize = tfile.RasterXSize
            dst_ysize = tfile.RasterYSize
            dst_gt = tfile.GetGeoTransform()
            tfile = None
            src = None
            
            # Get extent coordinates from reprojected raster
            xmin = dst_gt[0]
            ymax = dst_gt[3]
            xmax = xmin + dst_gt[1]*dst_xsize
            ymin = ymax + dst_gt[5]*dst_ysize               
                
            # Get columns and rows of new 30m dataset
            newxres = int((xmax - xmin) / xint) + 10
            newyres = int((ymax - ymin) / yint) + 10
             
            # Adjust xmin and ymax to align with the snap grid
            #newymax = (snap_y-(math.ceil(abs((ymax-snap_y)/30))*30))+30*5
            #newxmin = ((math.ceil(abs((xmin-snap_x)/30))*30)+snap_x)-30*5
            if ymax < snap_y:
                newymax = (snap_y-(math.ceil(abs((ymax-snap_y)/yint))*yint))+yint*5
            else:
                newymax = (snap_y+(math.ceil(abs((ymax-snap_y)/yint))*yint))+yint*5
            if xmin > snap_x:
                newxmin = (snap_x+(math.ceil(abs((xmin-snap_x)/xint))*xint))-xint*5
            else:
                newxmin = (snap_x-(math.ceil(abs((xmin-snap_x)/xint))*xint))-xint*5
                      
            newymin = newymax - newyres*yint
            newxmax = newxmin + newxres*xint
            
            # Set target raster to hold output
            memory_driver = gdal.GetDriverByName('GTiff')
            target_ds = memory_driver.Create(dst_filename, newxres, newyres, 1, gdal.GDT_Int16, options=['COMPRESS=LZW'])
        
            # Set the ROI image's projection and extent
            target_ds.SetProjection(rprj_wkt)
            target_ds.SetGeoTransform((newxmin, xint, 0, newymax, 0, -yint)) # use adusted coordinates for snapping
            # Fill output band with zeros
            #b = target_ds.GetRasterBand(1)
            #b.Fill(0)
            #b.FlushCache()
            target_ds = None
            
            pylogfor = odfor + '/pylog_' + k + '_' + continent + '.txt'
            with open(pylogfor, "a") as text_file:
                text_file.write("Warping " + src_filename + "\n")
                text_file.write("geotransform is " + str(dst_gt) + "\n")
                text_file.write("original column numbers is " + str(dst_xsize) + "\n")
                text_file.write("original row numbers is " + str(dst_ysize) + "\n")
                text_file.write("ymax is " + str(ymax) + "\n")
                text_file.write("ymin is " + str(ymin) + "\n")
                text_file.write("xmax is " + str(xmax) + "\n")
                text_file.write("xmin is " + str(xmin) + "\n")
                text_file.write("newymax is " + str(newymax) + "\n")
                text_file.write("newymin is " + str(newymin) + "\n")
                text_file.write("newxmax is " + str(newxmax) + "\n")
                text_file.write("newxmin is " + str(newxmin) + "\n")
            # Warp
            #subprocess.call(["gdalwarp", "-co", "COMPRESS=LZW", "-multi","-wo", "NUM_THREADS=ALL_CPUS", "-t_srs", rprj_wkt, "-r", "near", src_filename, dst_filename])
            subprocess.call(["gdalwarp", "-multi","-wo", "NUM_THREADS=ALL_CPUS", "-t_srs", rprj_wkt, "-r", "near", src_filename, dst_filename])
            # Convert to 8bit unsigned
            #transout = os.path.dirname(dst_filename) + os.sep + os.path.basename(dst_filename).split(".")[0] + "_rnd.tif"
            #subprocess.call(["gdal_translate", "-ot", "Byte",dst_filename, transout])
