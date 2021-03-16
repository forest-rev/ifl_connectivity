# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 13:08:31 2019
This script calculates cost distance corridors between point pairs and estimates effective
distances using cost distance values. It is a second stage script that expands
certain dimensions of the corridor mapping extent to avoid edge effects for particular
point pairs.
@author: pj276
"""

#%%
# Imports
import sys, os, ogr, math, shutil
import subprocess, gdal
from shapely.geometry import box
#from shapely.geometry import mapping
#import fiona
from skimage import graph
import pandas
import numpy as np
sys.path.append('/home/pj276/projects/ifl_corridors/code/functions/')
import hf_utilities as hfu
sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
import param_file as p1

#%%
# Set inputs/outputs
#headdir = p1.headdir
# Project directory
odir = p1.odir
# Finfile directory
ffdir = p1.ffdir
#'/projects/above_gedi/pjantz/ifl_corridors_data'
# Edge matched forest cover tiles folder
emd = p1.emd
#'/gfc/gfc1_4_sa_ss_emt'
# Corridor directory
sdir = p1.sdir
#headdir + '/outputs/regular_points/rp10ks1/corridors'
# Input hifo file
inpoly = p1.inpoly
#headdir + '/hifo/hland_tropics_ssu_sa.shp'
# Input csv folder and base name
if os.path.exists(p1.icsv2):
    # Set icsv to the set of redo files
    icsv = p1.icsv2
else:
    icsv = p1.icsv
#headdir + '/outputs/regular_points/rp10ks1/ptnbs_10k_300k_all.csv'
# list of cost surface files
forfileindex = p1.forfileindex
#headdir + emd + '/cs2000neg1index.txt'
# Name of redo file
rdname = p1.rdname
#'redome.txt'
# vrt file prefix
vrtp = p1.vrtp
#'/cs2000_neg1_'
# output corridor prefix
ocp = p1.ocp
#'/corr_neg1_'
# Output stats prefixe
optilep = p1.optilep
#'/lcc_ptile_dists_2000_neg1_'
xthp = p1.xthp
# Snap coordinates
snap_y = p1.snap_y
snap_x = p1.snap_x
# Pixel sizes (usually either 30 or 90m)
xint = p1.xint
yint = p1.yint
# Interaction distance
ixdist = p1.ixdist

#%%
# Functions to calculate dispersal probabilities
# Assume median dispersal distances of 100km, 10km, and 1km
# This function gets the value of alpha for a desired MEDIAN dispersal distance
def getmedexp(x):
  val = 1*(np.log(2))/x
  return(val)

#km100alpha = getmedexp(100) # 660km corresponds to a probability of 1%
#km25alpha = getmedexp(25) # 165km corresponds to a probability of 1%, 100km is 6.25%
#km10alpha = getmedexp(10) # 66km corresponds to a probability of 1%, 100km is 0.09%
#km1alpha = getmedexp(1) # 6.6km corresponds to a probability of 1%

# Use ixdist to calculate alpha
kmalpha = getmedexp(ixdist)

# Exponential distance decay
# x is distance, a is alpha
def incfun(x,a):
  op = np.exp(-x*a)
  return(op)

#%%
# Set parameters
# Job array index
k = int(sys.argv[1]) - 1 # subtract for python zero indexing
kstart = 250*k
kend = kstart + 250

# Read source target patches from csv columns to lists
# Columns are gc1,ptid1,cx1,cy1,gc2,ptid2,cx2,cy2
dat = pandas.read_csv(icsv)
# Select out desired columns
dat = dat[['gc1','ptid1','gc2','ptid2','cx1','cy1','cx2','cy2','dist']]

# Get subset of 250 (making sure to adjust end range if we're near the end of the table)
if kend > dat.shape[0]-1:
    kend = dat.shape[0] -1
dat = dat.loc[kstart:kend,]

# Group by gc1 and gc2
g1 = dat.groupby(['gc1', 'gc2'])

# Map each corridor in each group
for gname, g in g1:
    print gname
    # Create folder for each group (otherwise file system bogs down)
    tdir = sdir + '/' + str(int(gname[0])) + '_' + str(int(gname[1])) + '_' + str(k)
    if os.path.exists(tdir):
        shutil.rmtree(tdir)
        os.mkdir(tdir)
    else:
        os.mkdir(tdir)
    g = g.reset_index()
    # Create empty list to which corridor stat data frames will be appended
    cstatslist = []
    # Create empty list to which data frames giving corridor status (finished or invalid) will be appended
    # This is to keep track of which corridors have been processed
    cstatuslist = []
    # Create empty list to which data frames giving finished corridor names will be appended
    # This is to build the level0 tiff mosaic
    cfinishlist = []
    for j in range(0,g.shape[0]-1):
        print j
        print 'Corridor ID is ' + str(int(gname[0])) + '_' + str(int(gname[1])) + '_' + str(k) + '_' + str(j)
        # Read row to variables
        index,gc1,ptid1,gc2,ptid2,cx1,cy1,cx2,cy2,dist = g.loc[j, :].values.tolist()
        # Convert to ints
        gc1 = int(gc1)
        ptid1 = int(ptid1)
        gc2 = int(gc2)
        ptid2 = int(ptid2)    
        #%%
        # Set while loop variable
        expander = 1
        # Set counter for each boundary
        tecount = 0 # increment count by 0.5
        lecount = 0 # increment count by 0.5
        becount = 0 # increment count by 0.5
        recount = 0 # increment count by 0.5
        # Set boundary indicators to zero
        topedge = 0
        leftedge = 0
        bottomedge = 0
        rightedge = 0
        while expander == 1:
            #%%
            # Get min and max of point pair coordinates
            xmin = min(cx1,cx2); xmax = max(cx1,cx2)
            ymin = min(cy1,cy2); ymax = max(cy1,cy2)
            #%%
            # Use min and max extents from patches to create a circle with radius 150% of the bounding box width
            # Use bounding box of circle to get buffered extent for the cost surface
            # Halfway coordinate (buffer this coordinate by half the distance between points plus 25%)
            pdxmid = (xmax-xmin)/2+xmin
            pdymid = (ymax-ymin)/2+ymin
            # Half distance between points plus 50% or 75% or 100%.
            pdhalf = (((xmax-xmin)**2 + (ymax-ymin)**2)**0.5)/2*1.5
            # Convert to wkt
            wkt = "POINT (" + str(pdxmid) + " " + str(pdymid) + ")"
            # Create geometry
            pt = ogr.CreateGeometryFromWkt(wkt)
            # Buffer it
            poly = pt.Buffer(pdhalf)
            # Get extent of buffer. Use this to subset the cost raster
            bext = poly.GetEnvelope()
            bext = list(bext) # Tuple to list
            # If there's a boundary issue, add to the offending boundary(ies)
            # This set up should keep adding to edge extents if they ever had a boundary
            # issue which keeps them from resetting to the original envelope.
            if topedge == 1:
                bext[3] = bext[3] + abs((((bext[3]-bext[2])*(1+tecount))-(bext[3]-bext[2]))/2)
                #bext[3] = bext[3]*(1+tecount)
                tecount = tecount + 0.5
            else:
                bext[3] = bext[3] + abs((((bext[3]-bext[2])*(1+tecount))-(bext[3]-bext[2]))/2)
                #bext[3] = bext[3]*(1+tecount)
            if leftedge == 1:
                bext[0] = bext[0] - abs(((bext[1]-bext[0])*(1+lecount)-(bext[1]-bext[0]))/2)
                #bext[0] = bext[0]*(1+lecount)
                lecount = lecount + 0.5
            else:
                bext[0] = bext[0] - abs(((bext[1]-bext[0])*(1+lecount)-(bext[1]-bext[0]))/2)
                #bext[0] = bext[0]*(1+lecount)
            if bottomedge == 1:
                bext[2] = bext[2] - abs((((bext[3]-bext[2])*(1+becount))-(bext[3]-bext[2]))/2)
                #bext[2] = bext[2]*(1+becount)
                becount = becount + 0.5
            else:
                bext[2] = bext[2] - abs((((bext[3]-bext[2])*(1+becount))-(bext[3]-bext[2]))/2)
                #bext[2] = bext[2]*(1+becount)
            if rightedge ==1:
                bext[1] = bext[1] + abs(((bext[1]-bext[0])*(1+recount)-(bext[1]-bext[0]))/2)
                #bext[1] = bext[1]*(1+recount)
                recount = recount + 0.5
            else:
                bext[1] = bext[1] + abs(((bext[1]-bext[0])*(1+recount)-(bext[1]-bext[0]))/2)
                #bext[1] = bext[1]*(1+recount)
            print 'te-' + str(bext[3]) + ' le-' + str(bext[0]) + ' be-' + str(bext[2]) + ' re-' + str(bext[1])
            # Snap to forest cover grid
            newxres = int((bext[1] - bext[0]) / xint) + 1
            newyres = int((bext[3] - bext[2]) / yint) + 1
            if bext[3] < snap_y:
                newymax = (snap_y-(math.ceil(abs((bext[3]-snap_y)/yint))*yint))+yint*1
            else:
                newymax = (snap_y+(math.ceil(abs((bext[3]-snap_y)/yint))*yint))+yint*1
            if bext[0] > snap_x:
                newxmin = (snap_x+(math.ceil(abs((bext[0]-snap_x)/xint))*xint))-xint*1
            else:
                newxmin = (snap_x-(math.ceil(abs((bext[0]-snap_x)/xint))*xint))-xint*1  
            newymin = newymax - (newyres*yint)
            newxmax = newxmin + (newxres*xint)
            bext[0] = newxmin
            bext[1] = newxmax
            bext[2] = newymin
            bext[3] = newymax
            # Save extent to a new Shapefile
            boundingbox = box(bext[0], bext[2], bext[1], bext[3])
            #%%
            # Create resistance vrt tile
            csrast = tdir + vrtp + str(gc1) + '_' + str(ptid1) + '_' + str(gc2) + '_' + str(ptid2) + '.vrt'
            subprocess.call(["gdalbuildvrt", "-te", str(bext[0]), str(bext[2]), str(bext[1]), str(bext[3]), "-input_file_list", forfileindex, csrast])
            cso = hfu.raster2array(csrast)
            cs = cso[0]
            cs[cs == 255] = 1000 # set no go areas to higher value
            # Give x and y coordinates. Returns row and column indices of coordinates
            locs1 = hfu.coord2pixelOffset(csrast,[cx1],[cy1])
            locs1 = [i[0] for i in locs1] # drop list level
            locs1 = [(locs1[0],locs1[1])] # convert to a list of lists...
            # Get distance from point in first layer
            # Create graph object
            go = graph.MCP_Geometric(cs, offsets=None, sampling=(xint,yint), fully_connected=True)
            cc1,tb = go.find_costs(locs1)
            del tb # delete traceback
            locs2 = hfu.coord2pixelOffset(csrast,[cx2],[cy2])
            locs2 = [i[0] for i in locs2] # drop list level
            locs2 = [(locs2[0],locs2[1])] # convert to list of lists...
            # Create graph object
            go = graph.MCP_Geometric(cs, offsets=None, sampling=(xint,yint), fully_connected=True)
            cc2,tb = go.find_costs(locs2) # this may take x,y order...
            del tb # delete traceback
            # Add cost distance surfaces together
            cct = cc1+cc2
            cct = cct         
            # Calculate probability of dispersal values
            cct = incfun(cct/1000.0,kmalpha)
            # Check for all zeros (using float 32 bit depth.
            if np.max(cct.astype('float32')) == 0:
                cq10 = 0
            else:
                # Get highest xth percentile paths
                cq10 = np.percentile(cct[cct > 0], 100-xthp)
            #%%
            # Get lcp raster for invalid pixel check
            lcpraster = np.where(cct >= np.max(cct), 1, 0)
            # Check if corridor crosses invalid pixels or if it is all zeros.
            # If not, calculate stats
            if np.max(lcpraster*cs) != 1000 and cq10 != 0:
                # Convert to binary based on percentile threshold. 1=corridor, 0=not corridor
                # Check edges for ones, indicating that the extent needs expanding
                xx = np.where(cct >= cq10, 1, 0)
                # Set edge indicators 
                # Top edge
                if 1 in xx[0,0:xx.shape[1]-1]:
                    topedge = 1
                else:
                    topedge = 0
                # Left edge
                if 1 in xx[0:xx.shape[0]-1,0]:
                    leftedge = 1
                else:
                    leftedge = 0
                # Bottom edge
                if 1 in xx[xx.shape[0]-1,0:xx.shape[1]-1]:
                    bottomedge = 1
                else:
                    bottomedge = 0
                # Right edge
                if 1 in xx[0:xx.shape[0]-1,xx.shape[1]-1]:
                    rightedge = 1
                else:
                    rightedge = 0
                # Write any boundary issue to file.
                if 1 in xx[0,0:xx.shape[1]-1] or 1 in xx[0:xx.shape[0]-1,0] or 1 in xx[xx.shape[0]-1,0:xx.shape[1]-1] or 1 in xx[0:xx.shape[0]-1,xx.shape[1]-1]:
                    with open(tdir + '/expand_' + str(gc1) + '_' + str(gc2) + '_' + str(k) + '.txt', 'w') as text_file:
                        text_file.write(str(j) + '_' + str(gc1) + '_' + str(ptid1) + '_' + str(gc2) + '_' + str(ptid2) + '\n')
                else:
                    # Highest probability path value
                    lcd = np.max(cct[cct > 0])
                    # Average prob in corridor
                    apd = np.mean(cct[cct >= cq10])
                    # Get probability of dispersal values only within percentile corridor
                    cct = np.where(cct >= cq10, cct, 0)
                    # Output file name
                    orast = tdir + ocp + str(gc1) + '_' + str(ptid1) + '_' + str(gc2) + '_' + str(ptid2) + '_temp.tif'
                    #%%
                    # Save corridor raster
                    hfu.array2raster(orast, csrast, gdal.GDT_Float32, 'GTiff', cct)
                    # Append file
                    cfinishlist.append(orast)    
                    # Create new empty df and populate columns
                    newdf = pandas.DataFrame(columns=['gc1','ptid1','gc2','ptid2','cx1','cy1','cx2','cy2','dist','pdxmid','pdymid','lcd','apd'])
                    newdf = newdf.append({
                            'gc1': gc1,
                            'ptid1': int(ptid1),
                            'gc2': gc2,
                            'ptid2': int(ptid2),
                            'cx1': cx1,
                            'cy1': cy1,
                            'cx2': cx2,
                            'cy2': cy2,
                            'dist': dist,
                            'pdxmid': pdxmid,
                            'pdymid': pdymid,
                            'lcd': lcd,
                            'apd': apd}, ignore_index=True)
                    # Set index
                    newdf.index.name = 'rid'
                    # Append to list
                    cstatslist.append(newdf)
                    # Create new empty df and populate columns
                    findf = pandas.DataFrame(columns=['gc1','ptid1','gc2','ptid2','corridor', 'status'])
                    findf = findf.append({
                            'gc1': gc1,
                            'ptid1': int(ptid1),
                            'gc2': gc2,
                            'ptid2': int(ptid2),
                            'corridor': orast, 
                            'status': 'finished'}, ignore_index=True)
                    # Set index
                    findf.index.name = 'rid'
                    # Append to list
                    cstatuslist.append(findf)
                    os.remove(csrast)
                    # If there is no boundary issue, set expander to 0 and break out of while loop
                    expander = 0
            else:
                # Clean up
                os.remove(csrast)
                # Output raster name as if it had been processed    
                orast = tdir + ocp + str(gc1) + '_' + str(ptid1) + '_' + str(gc2) + '_' + str(ptid2) + '_temp.tif'
                # Create new empty df and populate columns
                findf = pandas.DataFrame(columns=['gc1','ptid1','gc2','ptid2','corridor', 'status'])
                findf = findf.append({
                        'gc1': gc1,
                        'ptid1': int(ptid1),
                        'gc2': gc2,
                        'ptid2': int(ptid2),
                        'corridor': orast, 
                        'status': 'invalid'}, ignore_index=True)
                # Set index
                findf.index.name = 'rid'
                # Append to list
                cstatuslist.append(findf)
                # Exit while loop
                expander = 0

    # Put in check here to make sure at least one file was processed
    if len(cstatuslist) > 0:
        # Concatenate data frames of corridor status
        findfall = pandas.concat(cstatuslist, ignore_index=True)
        # Subset finished corridors
        findfss = findfall.loc[findfall['status'] == 'finished']
        if findfss.shape[0] > 0:
            # Group by gc1 and gc2
            datgb = findfss.groupby(['gc1', 'gc2'])
            # Iterate through groups
            for name, group in datgb:
                # Build vrt of finished corridors
                ovrt = tdir + '/temp_vrt_' + str(int(name[0])) + '_' + str(int(name[1])) + '_' + str(k) + '.vrt'
                # Build a vrt from a file list
                acmd = ["gdalbuildvrt", "-separate", ovrt]
                for myfile in group['corridor']:
                    acmd.append(myfile)
                subprocess.call(acmd)
                
                # Iterative sum vrt to stack group corridors
                newRasterfn = tdir + ocp + str(int(name[0])) + '_' + str(int(name[1])) + '_' + str(k) + '.tif'
                zz = hfu.rastersum(ovrt)
                hfu.array2raster(newRasterfn, ovrt, gdal.GDT_Float32, 'GTiff', zz)
                
                # Delete individual files
                for myfile in group['corridor']:
                    os.remove(myfile)
                # Delete vrt
                os.remove(ovrt)
    
            # Concatenate data frames of corridor stats
            rdf = pandas.concat(cstatslist, ignore_index=True)
            # Write to csv
            rdf.to_csv(odir + optilep + '/cstat_' + str(int(name[0])) + '_' + str(int(name[1])) + '_' + str(k) + '.csv', index=False)
        
        # Write corridor status to csv
        # Group by gc1 and gc2
        g2 = findfall.groupby(['gc1', 'gc2'])
        for name, group in g2:
            findf.to_csv(odir + ffdir + '/finstatus_' + str(int(name[0])) + '_' + str(int(name[1])) + '_' + str(k) + '.csv', index=False)
