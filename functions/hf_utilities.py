# -*- coding: utf-8 -*-
"""
Created on Wed May 10 14:48:28 2017

@author: pj276
"""

#---------------------------------------------------------------------------------------------
# This file holds functions that, among other things, are useful for reading
# and manipulating spatial data using gdal.
#---------------------------------------------------------------------------------------------

# Import libraries
import gdal, ogr, osr, math, os, sys, urllib, itertools, pickle
import geopandas as gpd
import fiona
#from skimage import graph
#from scipy import ndimage as ndi
from scipy.ndimage.filters import maximum_filter, generic_filter
from skimage import graph
from skimage.morphology import binary_erosion, binary_dilation
#import mahotas
import numpy as np
import shapefile as shp
from scipy.sparse import csr_matrix
from scipy.sparse import identity
from scipy.sparse import diags
import scipy.sparse.linalg as spl

# Define functions
# Save an object to disk as a pkl file
def save_object(obj, filename):
    ''' # sample usage
    save_object(company1, 'company1.pkl')
    '''
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
'''
Provide raster name as a string. Provide clip as a list if using.
Clip order is xmin, ymin, xmax, ymax
Geotransform elements are originX, pixelWidth, 0, originY, 0, pixelHeight
pixelHeight is generally negative so multiply by -1 when calculating ymin
Clip takes the intersection of the extents
'''
def raster2array(rasterfn, clip=False):
    if clip != False:
        raster = gdal.Open(rasterfn)
        geotransform = raster.GetGeoTransform()
        clip[0] = max(clip[0], geotransform[0])
        clip[1] = max(clip[1], geotransform[3]-np.float((geotransform[5]*-1*raster.RasterYSize)))
        clip[2] = min(clip[2], geotransform[0]+np.float((geotransform[1]*raster.RasterXSize)))
        clip[3] = min(clip[3], geotransform[3])
        inv_geotransform = gdal.InvGeoTransform(geotransform)
        _x0, _y0 = gdal.ApplyGeoTransform(inv_geotransform, clip[0], clip[1])
        _x1, _y1 = gdal.ApplyGeoTransform(inv_geotransform, clip[2], clip[3])
        x0, y0 = min(_x0, _x1), min(_y0, _y1)
        x1, y1 = max(_x0, _x1), max(_y0, _y1)
        band = raster.GetRasterBand(1)
        array = band.ReadAsArray(int(np.round(x0)), int(np.round(y0)), int(np.round(x1-x0)), int(np.round(y1-y0)))
        try:
            rysize = array.shape[0]
            rxsize = array.shape[1]
        except AttributeError:
            print('Empty array, no overlap between subset and original raster')
            return(False)
        proj = raster.GetProjection()
        geotransform = (clip[0], geotransform[1], geotransform[2], clip[3], geotransform[4], geotransform[5])
        return array,geotransform,rxsize,rysize,proj
    else:
        raster = gdal.Open(rasterfn)
        band = raster.GetRasterBand(1)
        array = band.ReadAsArray()
        geotransform = raster.GetGeoTransform()
        rxsize = raster.RasterXSize
        rysize = raster.RasterYSize
        proj = raster.GetProjection()
        return array,geotransform,rxsize,rysize,proj

def raster2array3d(rasterfn):
    raster = gdal.Open(rasterfn)
    array = raster.ReadAsArray()
    geotransform = raster.GetGeoTransform()
    rxsize = raster.RasterXSize
    rysize = raster.RasterYSize
    proj = raster.GetProjection()
    return array,geotransform,rxsize,rysize,proj

def rastersum(rasterfn):
    '''This function takes a 3d vrt and iteratively sums the bands. It outputs
    an array.
    '''
    raster = gdal.Open(rasterfn)
    rsum = raster.GetRasterBand(1)
    rsum = rsum.ReadAsArray().astype('int32')
    rcount = raster.RasterCount
    if rcount > 1:
        for band in range(1,rcount):
            band += 1
            rb = raster.GetRasterBand(band)
            if rb is None:
                continue
            ra = rb.ReadAsArray()
            rsum = rsum + ra
    return rsum

def rastersumcfx(rasterfn, distlist):
    '''This function takes a 3d vrt and iteratively sums the bands. It outputs
    an array. It also applies a multiplier from list elements to 
    each band. It outputs two rasters. One is the sum of the binary rasters
    and the other is the sum of rasters with multiplier applied
    '''
    raster = gdal.Open(rasterfn)
    rsum = raster.GetRasterBand(1)
    rsum1 = rsum.ReadAsArray().astype('uint32')
    rsum2 = rsum1*distlist[0]
    rcount = raster.RasterCount
    if rcount > 1:
        for band in range(2,rcount):
            #band += 1
            rb = raster.GetRasterBand(band)
            if rb is None:
                continue
            ra = rb.ReadAsArray()
            rsum1 = rsum1 + ra
            acd = distlist[band-1]
            ra = ra*acd
            rsum2 = rsum2 + ra
    return [rsum1, rsum2]

def coord2pixelOffset(rasterfn,x,y):
    ''' This function converts coordinates to matrix indices. It takes
    a raster file name and a list of x coordinates and a list of y coordinates.
    It returns a tuple of lists of row and column indices.
    '''
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0] # I think this is upper left corner
    originY = geotransform[3] # I think this is upper left corner
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    xOffset = [int((i - originX)/pixelWidth) for i in x]
    yOffset = [int((i - originY)/pixelHeight) for i in y]
    return yOffset,xOffset

# Need to double check that this gets the right pixel
# This returns the coordinate of the pixel center
def pixel2coordOffsetcc(rasterfn,x,y):
    ''' This function converts matrix indices to coordinates. It takes
    a raster file name and a list of x indices and a list of y indices.
    It returns a list of x and a list of y coordinates.
    '''
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    xco = [originX+(i*pixelWidth) + (0.5*pixelWidth) for i in x]
    yco = [originY+(i*pixelHeight) + (0.5*pixelHeight) for i in y]
    return xco,yco

def array2raster(newRasterfn,rasterfn,ptype,gdriver,array,geotrans=False,rasterproj=False):
    ''' This function writes an array to tiff. It uses an existing tiff
    as a template or a geotransform and wkt projection. 'array' is given as an object. 
    newRasterfn=new raster file name as a string.
    rasterfn=name of existing raster to use as a template as a string.
    ptype=pixel type of output as gdal object (e.g. gdal.GDT_Byte or gdal.GDT_Float32)
    gdriver=string giving gdal driver type (e.g.'GTiff', 'VRT')
    array=numpy array
    geotrans=list
    rasterproj=well known text
    '''
    if geotrans != False:
        geotransform = geotrans
        originX = geotransform[0]
        originY = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]
        cols = array.shape[1]
        rows = array.shape[0]
        driver = gdal.GetDriverByName(gdriver)
        outRaster = driver.Create(newRasterfn, cols, rows, 1, ptype, ['COMPRESS=LZW'])
        outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
        outband = outRaster.GetRasterBand(1)
        outband.WriteArray(array)
        outRasterSRS = osr.SpatialReference()
        outRasterSRS.ImportFromWkt(rasterproj)
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        outband.FlushCache()
    else:
        raster = gdal.Open(rasterfn)
        geotransform = raster.GetGeoTransform()
        originX = geotransform[0]
        originY = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]
        cols = raster.RasterXSize
        rows = raster.RasterYSize
        driver = gdal.GetDriverByName(gdriver)
        outRaster = driver.Create(newRasterfn, cols, rows, 1, ptype, ['COMPRESS=LZW'])
        outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
        outband = outRaster.GetRasterBand(1)
        outband.WriteArray(array)
        outRasterSRS = osr.SpatialReference()
        outRasterSRS.ImportFromWkt(raster.GetProjectionRef())
        outRaster.SetProjection(outRasterSRS.ExportToWkt())
        outband.FlushCache()


def array2rasterclip(newRasterfn,rasterfn,ptype,gdriver,array,xoff,yoff):
    ''' This function writes an array to tiff. It uses an existing tiff
    as a template but clips to the extent of the array. Useful for
    saving subsets.
    'array' is given as an object. 
    newRasterfn=new raster file name as a string.
    rasterfn=name of existing raster to use as a template as a string.
    ptype=pixel type of output as gdal object (e.g. gdal.GDT_Byte or gdal.GDT_Float32)
    gdriver=string giving gdal driver type (e.g.'GTiff', 'VRT')
    array=numpy array
    xoff=xoffset in meters, needs to be multiple of pixel size
    yoff=yoffset in meters, needs to be multiple of pixel size
    Give offsets as negative if new array is smaller than original raster
    Give offsets as positive if new array is larger than original raster
    '''
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]-xoff
    originY = geotransform[3]+yoff
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    cols = array.shape[1]
    rows = array.shape[0]
    driver = gdal.GetDriverByName(gdriver)
    outRaster = driver.Create(newRasterfn, cols, rows, 1, ptype, ['COMPRESS=LZW'])
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(raster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

# Function to return raster corner coordinates from geotransform and pixel size
def raster_corners(srcdataset):
    src = gdal.Open(srcdataset)
    gt = src.GetGeoTransform()
    xsize = src.RasterXSize
    ysize = src.RasterYSize
    xmin = gt[0]
    ymin = gt[3] + ysize*gt[5]
    xmax = gt[0] + xsize*gt[1]
    ymax = gt[3]
    return([xmin,ymin,xmax,ymax])
    
# Function to return neighboring regions from a footprint (aka structuring element)
# Used as input for list_neighbors function
def region_list(values, elist):
    values = set(values.astype(int)) # Set values to integer
    values = [x for x in values if x != -1] # Remove background value
    if len(values) > 1:
        if tuple(values) not in elist and tuple(values)[::-1] not in elist:
            elist.append(tuple(values)) # Add adjacent regions to list
    return 0.0 # Return zero for filter array output

# Function to return a list of neighbors. This function takes 'region_list'
# as the filter function input. The footprint (aka structuring element) is
# defined in the function
def list_neighbors(thiessen_poly,filter_fun):
    fp = np.array([[1, 1]])
    hlist=[]
    tarray = generic_filter(thiessen_poly,
                           function=filter_fun,
                           footprint=fp,
                           mode='nearest',
                           extra_arguments=(hlist,))
    del tarray
    # Get neighboring regions (vertical neighbors)
    vlist = []
    fp = np.array([[1],
                    [1]])
    tarray = generic_filter(thiessen_poly,
                           function=filter_fun,
                           footprint=fp,
                           mode='nearest',
                           extra_arguments=(vlist,))
    del tarray
    # Combine lists of neighbors, using hlist as the base list
    addregions = [i for i in vlist if i not in hlist and i[::-1] not in hlist]
    for i in addregions:
        hlist.append(i)
    return(hlist)

# Function to get least cost path from each cell to the nearest origin.
# Do this by iterating through the cumulative cost array and calling traceback each time.
# The result is a raster of thiessen polygon zones
def tpfunction(ccost,mcp_obj):
    tp = np.zeros(ccost.shape)
    it = np.nditer(ccost, flags=['multi_index'])
    while not it.finished:
        print it.multi_index
        # Get row and column of current cell
        arow = it.multi_index[0] # row
        acol = it.multi_index[1] # column
        # If cumulative value is infinite, set to -1
        if np.isinf(ccost[arow,acol]):
            tp[arow,acol] = -1
            it.iternext()
        # Otherwise, get index of closest origin
        else:
            it[0] # array value
            apath = mcp_obj.traceback([arow,acol])
            co = apath[0] # row,col of closest origin
            coi = (co[0]-1)*(co[1]-1) # index of closest origin (1-indexed, row first, then column)
            # Update zeros array with coi
            tp[arow,acol] = coi
            it.iternext()
    # Convert to integer
    tp = tp.astype(int)
    return(tp)

# This function finds a least cost path through an array between supplied indices
# Modified from - https://gis.stackexchange.com/questions/28583/gdal-perform-simple-least-cost-path-analysis/78699#78699
# startIndex is supplied as a list with two elements. Same with stopIndex.
# Output is a list, first element is a binary array, second is a list of indices
def lcpfunction(costSurfaceArray,startIndex,stopIndex):   
    # coordinates to array index
    startIndexY,startIndexX = startIndex[0],startIndex[1]
    stopIndexY,stopIndexX = stopIndex[0],stopIndex[1]
    # create path
    indices, weight = graph.route_through_array(costSurfaceArray, (startIndexY,startIndexX), (stopIndexY,stopIndexX),geometric=True,fully_connected=True)
    indices = np.array(indices).T
    path = np.zeros_like(costSurfaceArray)
    path[indices[0], indices[1]] = 1
    return path,indices

# Function to get great circle distance
def distance(origin, destination):
    lat1, lon1 = origin
    lat2, lon2 = destination
    radius = 6371 # km
    dlat = math.radians(lat2-lat1)
    dlon = math.radians(lon2-lon1)
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c
    return d

# This function coverts coordinates to point Shapefile
# It takes a list of tuples with x,y ordering (coordlist)
# Outfile (ofname) should include .shp extension
# User supplies a projection code (epsgcode) as a string e.g. "4326"
# or as a SpatialReference class from osgeo.osr
def coordinates2point(coordlist,ofname,sr):
    # Instantiate shapefile writer
    w = shp.Writer(shp.POINT)
    # Add field
    w.field('id','N')
    for i,j in enumerate(coordlist):
        # Input data
        w.point(float(j[0]),float(j[1])) # -124.4577,48.0135
        w.record(i)
    w.save(ofname)
    if isinstance(sr,str):
        # Add projection information (from https://glenbambrick.com/2015/08/09/prj/)
        wkt = urllib.urlopen("http://spatialreference.org/ref/epsg/{0}/prettywkt/".format(sr))
        # remove spaces between characters
        remove_spaces = wkt.read().replace(" ","")
        # place all the text on one line
        output = remove_spaces.replace("\n", "")
    elif str(type(sr)) == "<class 'osgeo.osr.SpatialReference'>":
        sr.MorphToESRI()
        output = sr.ExportToWkt()
    # Replace shp extension with prj extension
    prjname = ofname.replace(".shp",".prj")
    prj = open(prjname, "w")
    prj.write(output)
    prj.close()

# This function takes a polygon layer as an input and converts it to
# a line.
def poly2line(input_poly,output_line):
    #source_ds = ogr.Open(input_poly)
    #source_layer = source_ds.GetLayer()
    source_layer = input_poly
    # get spatial ref from input poly layer
    sr = source_layer.GetSpatialRef()
    # polygon2geometryCollection
    geomcol =  ogr.Geometry(ogr.wkbGeometryCollection)
    for feat in source_layer:
        geom = feat.GetGeometryRef()
        ring = geom.GetGeometryRef(0)
        geomcol.AddGeometry(ring)
    # geometryCollection2shp
    shpDriver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(output_line):
            shpDriver.DeleteDataSource(output_line)
    outDataSource = shpDriver.CreateDataSource(output_line)
    outLayer = outDataSource.CreateLayer(output_line, geom_type=ogr.wkbMultiLineString, srs=sr)
    featureDefn = outLayer.GetLayerDefn()
    outFeature = ogr.Feature(featureDefn)
    outFeature.SetGeometry(geomcol)
    outLayer.CreateFeature(outFeature)
    outFeature = None

# This function takes a line as input and converts it to a raster.
# xint is x pixel size. yint is y pixel size. burnfield is the
# a string of a field name to use for raster values.
# Options for rasterize layer needs to be a list e.g. ['ATTRIBUTE=' + burnfield,'ALL_TOUCHED=TRUE']
def line2raster(inline,outraster,xint,yint,snap_y,snap_x,burnconstant=1000,options=None):
    # Convert line to raster
    # Define pixel_size and NoData value of new raster
    #pixel_size = 0.000277777777777778
    # Filename of input OGR file
    vector_fn = inline
    #vector_fn = "C:/Users/pj276/Projects/ifl_corridors/geodata/hfl_2.shp"
    # Open the data source and read in the extent
    source_ds = ogr.Open(vector_fn)
    source_layer = source_ds.GetLayerByIndex(0)
    x_min, x_max, y_min, y_max = source_layer.GetExtent()
    #geometry = source_layer.GetGeomType()
    #geometry_name = ogr.GeometryTypeToName(geometry)
    #print("The layer's geometry is: {geom}\n".format(geom=geometry_name))
    # Get projection info
    sr = source_layer.GetSpatialRef()
    sr = sr.ExportToWkt() # convert to well known text format
    # Create the destination data source. Buffer by 10 pixels on each side.
    x_res = int((x_max - x_min) / xint) + 10
    y_res = int((y_max - y_min) / yint) + 10
    # Calculate a new upper left coordinate to align with snap grid
    # Add half of buffer to top and left side of grid
    cpy = (snap_y-(math.ceil(abs((y_max-snap_y)/30))*30))+30*5
    cpx = ((math.ceil(abs((x_min-snap_x)/30))*30)+snap_x)-30*5
          
    # Adjust xmin and ymax to align with the snap grid
    #newymax = (snap_y-(math.ceil(abs((ymax-snap_y)/30))*30))+30*5
    #newxmin = ((math.ceil(abs((xmin-snap_x)/30))*30)+snap_x)-30*5
    if y_max < snap_y:
        cpy = (snap_y-(math.ceil(abs((y_max-snap_y)/yint))*yint))+yint*5
    else:
        cpy = (snap_y+(math.ceil(abs((y_max-snap_y)/yint))*yint))+yint*5
    if x_min > snap_x:
        cpx = (snap_x+(math.ceil(abs((x_min-snap_x)/xint))*xint))-xint*5
    else:
        cpx = (snap_x-(math.ceil(abs((x_min-snap_x)/xint))*xint))-xint*5  

    # Set target raster to hold output
    memory_driver = gdal.GetDriverByName('GTiff')
    target_ds = memory_driver.Create(outraster, x_res, y_res, 1, gdal.GDT_Int16)
    #target_ds = memory_driver.Create('C:/Users/pj276/Projects/ifl_corridors/geodata/hfb_2.tif', x_res, y_res, 1, gdal.GDT_Byte)
    # Set the ROI image's projection and extent
    target_ds.SetProjection(sr)
    target_ds.SetGeoTransform((cpx, xint, 0, cpy, 0, -xint)) # use adusted coordinates for snapping
    # Fill output band with zeros
    b = target_ds.GetRasterBand(1)
    b.Fill(0)
    # Rasterize
    status = gdal.RasterizeLayer(target_ds,  # output to our new dataset
                                 [1],  # output to our new dataset's first band
                                 source_layer,  # rasterize this layer
                                 None, None,  # don't worry about transformations since we're in same projection
                                 [burnconstant], # Code pixels with this value
                                 options # example ['ATTRIBUTE=' + burnfield,'ALL_TOUCHED=TRUE']
                                 #[100], # Code line pixels as 100
                                 )
    if status != 0:
        print("Conversion failed")
    else:
        print("Success")
    # Close dataset
    target_ds = None

# This function converts raster array coordinates from their original
# projection to one specified by a Basemap object
def convertXYraster(xy_source, inproj, outproj):
    # Get dimensions
    shape = xy_source[0,:,:].shape
    size = xy_source[0,:,:].size
    # the ct object takes and returns pairs of x,y, not 2d grids
    # so the the grid needs to be reshaped (flattened) and back.
    ct = osr.CoordinateTransformation(inproj, outproj)
    xy_target = np.array(ct.TransformPoints(xy_source.reshape(2, size).T))
    xx = xy_target[:,0].reshape(shape)
    yy = xy_target[:,1].reshape(shape)
    return xx, yy

# This function converts coordinates from source to target projections
# inproj is given as an osr spatial reference object
# outproj is given as well known text
#import osr, gdal    
#Get projection from raster and convert to wkt
#rwkt = gdal.Open(dfrname)
#rwkt = rwkt.GetProjectionRef()
#Get projection from epsg number
#swkt = osr.SpatialReference()
#swkt.ImportFromEPSG(4326)
def convertXYcoords(xco, yco, inproj, outproj):
    outproj=osr.SpatialReference(wkt=outproj)
    ct = osr.CoordinateTransformation(inproj, outproj)
    x,y,z = ct.TransformPoint(xco,yco)
    return x,y

def reproPoints(pt_shape, outproepsg, fsave='n'):
    """ Reprojects coordinates of a point shapefile.
    Takes a full file path to a point shapefile, an EPSG string (e.g. '4326'), and 'y' or 'n'
    indicating whether to save the output. Defaults to 'n'. In either case, it returns a
    shapefile writer object that can be used to grab the transformed points. E.g.,if output
    is assigned to swo, [i.points for i in swo.shapes()] returns a list of transformed 
    coordinates of the shape objects. Modified from
    https://glenbambrick.com/2016/01/24/reproject-shapefile/
    """
    try:
        #------------------------------------
        # Get original shapefile projection
        prjf = pt_shape.replace('.shp','.prj') # Get projection file name
        f = open(prjf) # Open projection file
        fr = f.read() # Read to variable
        srs = osr.SpatialReference(wkt=fr) # Convert to osr spatial reference
        f = None # Close file
        #------------------------------------
        # Create spatial reference for new shapefile (from https://glenbambrick.com/2015/08/09/prj/)
        wkt = urllib.urlopen("http://spatialreference.org/ref/epsg/{0}/prettywkt/".format(outproepsg))
        # remove spaces between charachters
        remove_spaces = wkt.read().replace(" ","")
        # place all the text on one line
        output = remove_spaces.replace("\n", "")
        srsnew = osr.SpatialReference(wkt=output) # Spatial reference for new shapefile
        # Create coordinate transform object
        ct = osr.CoordinateTransformation(srs, srsnew)
        #------------------------------------
        # Read in shapefile
        sf = shp.Reader(pt_shape)
        newshp = shp.Writer(shp.POINT) # Create object for new shapefile
        fields = sf.fields # Get shapefile fields
        #newfields = newshp.fields # Create object for new shapefile fields
        #------------------------------------
        # Transfer field structure to new shapefile
        for name in fields:
            if type(name) == "tuple":
                continue
            else:
                args = name
                newshp.field(*args)
        #------------------------------------
        # Transfer records to new shapefile fields
        records = sf.records()
        for row in records:
            args = row
            newshp.record(*args)
        #------------------------------------
        geom = sf.shapes() # Access original shapefile geometry
        for feature in geom:
            coords = feature.points[0]
            x, y = coords[0], coords[1]
            # transform the coordinates
            new_x,new_y,new_z = ct.TransformPoint(x,y)
            newshp.point(new_x,new_y)
        if fsave == 'y':
            oname = pt_shape.replace('.shp','')
            oname = oname + '_repro.shp'
            newshp.save(oname)
            # Write projection file
            oprj = oname.replace('.shp','.prj')
            prj = open(oprj, "w")
            prj.write(output)
            prj.close()
        return newshp
    except IOError:
        print 'Projection file does not appear to exist. Aborting.'

def getPointCoords(pt_shape):
    sf = shp.Reader(pt_shape) # Read in shapefile
    shapes = sf.shapes() # Create shapes object
    clist = [i.points for i in shapes] # Get point coordinates
    clist = list(itertools.chain.from_iterable(clist)) # Remove one list level
    xco = [i[0] for i in clist] # Get x coordinate
    yco = [i[1] for i in clist] # Get y coordinate
    return xco,yco

def getPointRecords(pt_shape):
    sf = shp.Reader(pt_shape) # Read in shapefile
    recs = sf.shapeRecords() # Create records object
    rlist = [i.record for i in recs] # Get records
    rlist = list(itertools.chain.from_iterable(rlist)) # Remove one list level
    return rlist

def findLocMax(in_array,fp):
    """
    Takes an array and detect the peaks usingthe local maximum filter.
    Returns a boolean mask of the peaks (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    # This is from Ivan's answer to a question on Stack Exchange
    # https://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array
    # I made some minor changes.    
    #apply the local maximum filter; all pixel of maximal value 
    #in their neighborhood are set to 1
    # I think the algorithm takes a derivative of the surface
    local_max = maximum_filter(in_array, footprint=fp)==in_array
    #local_max is a mask that contains the peaks we are 
    #looking for, but also the background.
    #In order to isolate the peaks we must remove the background from the mask.
    #we create the mask of the background
    background = (in_array==0)
    #a little technicality: we must erode the background in order to 
    #successfully subtract it form local_max, otherwise a line will 
    #appear along the background border (artifact of the local maximum filter)
    eroded_background = binary_erosion(background, structure=fp, border_value=1)
    #we obtain the final mask, containing only peaks, 
    #by removing the background from the local_max mask (xor operation)
    dps = local_max ^ eroded_background
    return(dps)

def reproLines(line_shape, sr, fsave='n'):
    """ Reprojects coordinates of a polyline shapefile.
    Takes a full file path to a polyline shapefile, an EPSG string (e.g. '4326') or osr spatial 
    reference object, and 'y' or 'n' indicating whether to save the output. Defaults to 'n'.
    In either case, it returns a shapefile writer object that can be used to grab the transformed data.
    Modified from https://glenbambrick.com/2016/01/24/reproject-shapefile/
    """
    try:
        #------------------------------------
        # Get original shapefile projection
        prjf = line_shape.replace('.shp','.prj') # Get projection file name
        f = open(prjf) # Open projection file
        fr = f.read() # Read to variable
        srs = osr.SpatialReference(wkt=fr) # Convert to osr spatial reference
        f = None # Close file
        #------------------------------------
        # Get new projection
        if isinstance(sr,str):
            # Add projection information (from https://glenbambrick.com/2015/08/09/prj/)
            wkt = urllib.urlopen("http://spatialreference.org/ref/epsg/{0}/prettywkt/".format(sr))
            # remove spaces between characters
            remove_spaces = wkt.read().replace(" ","")
            # place all the text on one line
            output = remove_spaces.replace("\n", "")
            srsnew = osr.SpatialReference(wkt=output) # Spatial reference for new shapefile
        elif str(type(sr)) == "<class 'osgeo.osr.SpatialReference'>":
            sr.MorphToESRI()
            output = sr.ExportToWkt()
            srsnew = sr # Spatial reference for new shapefile
        # Create coordinate transform object
        ct = osr.CoordinateTransformation(srs, srsnew)
        #------------------------------------
        # Read in shapefile
        sf = shp.Reader(line_shape)
        newshp = shp.Writer(shp.POLYLINE) # Create object for new shapefile
        fields = sf.fields # Get shapefile fields
        #------------------------------------
        # Transfer field structure to new shapefile
        for name in fields:
            if type(name) == "tuple":
                continue
            else:
                args = name
                newshp.field(*args)
        #------------------------------------
        # Transfer records to new shapefile fields
        records = sf.records()
        for row in records:
            args = row
            newshp.record(*args)
        #------------------------------------
        geom = sf.shapes() # Access original shapefile geometry
        for feature in geom:
            # if there is only one part
            if len(feature.parts) == 1:
                # create empty list to store all the coordinates
                poly_list = []
                # get each coord that makes up the polygon
                for coords in feature.points:
                    x, y = coords[0], coords[1]
                    # transform the coord
                    new_x,new_y,new_z = ct.TransformPoint(x,y)
                    # put the coord into a tuple structure
                    poly_coord = (float(new_x), float(new_y))
                    # append the coords to the polygon list
                    poly_list.append(poly_coord)
                # add the geometry to the shapefile.
                newshp.line(parts=[poly_list])
            else:
                print "this is a multi-part polyline"
                # append the total amount of points to the end of the parts list
                feature.parts.append(len(feature.points))
                # enpty list to store all the parts that make up the complete feature
                poly_list = []
                # keep track of the part being added
                parts_counter = 0
                # while the parts_counter is less than the amount of parts
                while parts_counter < len(feature.parts) - 1:
                    # keep track of the amount of points added to the feature
                    coord_count = feature.parts[parts_counter]
                    # number of points in each part
                    no_of_points = abs(feature.parts[parts_counter] - feature.parts[parts_counter + 1])
                    # create list to hold individual parts - these get added to poly_list[]
                    part_list = []
                    # cut off point for each part
                    end_point = coord_count + no_of_points
                    # loop through each part
                    while coord_count < end_point:
                        for coords in feature.points[coord_count:end_point]:
                            x, y = coords[0], coords[1]
                            # transform the coord
                            new_x,new_y,new_z = ct.TransformPoint(x,y)
                            # put the coord into a tuple structure
                            poly_coord = (float(new_x), float(new_y))
                            # append the coords to the part list
                            part_list.append(poly_coord)
                            coord_count = coord_count + 1
                    # append the part to the poly_list
                    poly_list.append(part_list)
                    parts_counter = parts_counter + 1
                # add the geometry to to new file
                newshp.line(parts=poly_list) # note this is slightly different than the single part case
        if fsave == 'y':
            oname = line_shape.replace('.shp','')
            oname = oname + '_repro.shp'
            newshp.save(oname)
            # Write projection file
            oprj = oname.replace('.shp','.prj')
            prj = open(oprj, "w")
            prj.write(output)
            prj.close()
        return newshp
    except IOError:
        print 'Projection file does not appear to exist. Aborting.'

# Check feature intersection using ogr
# Mostly from https://github.com/csterling/Morps/blob/9d16ca18a3fc952b1d7436322d7704bcfa04d9ff/morps.py
# f1 is the fishnet shapefile
# f2 is the tile index
def intersect(f1,f2,fid1=0,fid2=0):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    file1 = driver.Open(f1,0)
    layer1 = file1.GetLayer()
    feat1 = layer1.GetFeature(fid1)
    geom1 = feat1.GetGeometryRef()
    file2 = driver.Open(f2,0)
    layer2 = file2.GetLayer()
    feat2 = layer2.GetFeature(fid2)
    geom2 = feat2.GetGeometryRef()
    if geom1.Intersects(geom2) == 1:
        return feat2.GetFieldAsString(0)
    else:
        return None

# This function calculates the area of polygon intersection
def intersectionArea(f1,f2,fid1=0,fid2=0):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    file1 = driver.Open(f1,0)
    layer1 = file1.GetLayer()
    feat1 = layer1.GetFeature(fid1)
    geom1 = feat1.GetGeometryRef()
    file2 = driver.Open(f2,0)
    layer2 = file2.GetLayer()
    feat2 = layer2.GetFeature(fid2)
    geom2 = feat2.GetGeometryRef()
    if geom1.Intersects(geom2) == 1:
        intsect = geom1.Intersection(geom2)
        reachableArea = intsect.GetArea()
        return reachableArea
    else:
        return None

# Get a feature count from a shapefile
def getfc(f1):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    file1 = driver.Open(f1,0)
    layer1 = file1.GetLayer()
    fc = layer1.GetFeatureCount()
    return fc

# Create a regular fishnet polygon
# Edit so that the function takes upper left corner as a starting coordinate.
# It takes the number of rows and columns as well.
# Adapted from https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html#create-fishnet-grid
# fishnet(fnetfn, "-10464833.986991045996547", "-4185773.986991046", "-3320278.3702203743", "3320251.629779625684023", "999990", "999990")
# I believe outputShapefn needs a file extension e.g. .shp
# Arguments all given as strings and in gridHeight and gridWidth in map units
def fishnet(outputShapefn,xmin,ymax,rows,cols,gridHeight,gridWidth):
    # convert sys.argv to numeric
    xmin = float(xmin)
    #xmax = float(xmax)
    #ymin = float(ymin)
    ymax = float(ymax)
    cols = int(cols)
    rows = int(rows)
    gridWidth = float(gridWidth)
    gridHeight = float(gridHeight)
    # Calculate xmax
    #xmax = (gridWidth*cols)+xmin
    # Calculate ymin
    #ymin = (gridHeight*rows)-ymax
    # get rows
    #rows = math.ceil((ymax-ymin)/gridHeight)
    # get columns
    #cols = math.ceil((xmax-xmin)/gridWidth)
    # start grid cell envelope
    ringXleftOrigin = xmin
    ringXrightOrigin = xmin + gridWidth
    ringYtopOrigin = ymax
    ringYbottomOrigin = ymax-gridHeight
    # create output file
    outDriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(outputShapefn):
        os.remove(outputShapefn)
    outDataSource = outDriver.CreateDataSource(outputShapefn)
    outLayer = outDataSource.CreateLayer(outputShapefn,geom_type=ogr.wkbPolygon )
    featureDefn = outLayer.GetLayerDefn()
    # create grid cells
    countcols = 0
    while countcols < cols:
        countcols += 1
        # reset envelope for rows
        ringYtop = ringYtopOrigin
        ringYbottom =ringYbottomOrigin
        countrows = 0
        while countrows < rows:
            countrows += 1
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYtop)
            ring.AddPoint(ringXrightOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYbottom)
            ring.AddPoint(ringXleftOrigin, ringYtop)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)
            # add new geom to layer
            outFeature = ogr.Feature(featureDefn)
            outFeature.SetGeometry(poly)
            outLayer.CreateFeature(outFeature)
            outFeature = None
            # new envelope for next poly
            ringYtop = ringYtop - gridHeight
            ringYbottom = ringYbottom - gridHeight
        # new envelope for next poly
        ringXleftOrigin = ringXleftOrigin + gridWidth
        ringXrightOrigin = ringXrightOrigin + gridWidth
    # Save and close DataSources
    outDataSource = None

#print "[ ERROR ] you must supply seven arguments: output-shapefile-name.shp xmin xmax ymin ymax gridHeight gridWidth"

# This function takes a multiband input vrt and makes a maximum or sum mosaic
# of the vrt within the extent of a polygon in an input fishnet shapefile.
# It also needs an index for selecting a feature from the fishnet.
# Mosfun should be either 'sum' or 'max'
def tile_mosaic(input_zone_polygon, input_value_raster, newRasterfn, ind, mosfun):
    # Open data
    raster = gdal.Open(input_value_raster)
    shp = ogr.Open(input_zone_polygon)
    lyr = shp.GetLayer()
    # Get raster georeference info
    transform = raster.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = transform[5]
    # Get vector geometry
    feat = lyr.GetFeature(ind)
    # Get extent of feat
    geom = feat.GetGeometryRef()
    if (geom.GetGeometryName() == 'MULTIPOLYGON'):
        count = 0
        pointsX = []; pointsY = []
        for polygon in geom:
            geomInner = geom.GetGeometryRef(count)
            ring = geomInner.GetGeometryRef(0)
            numpoints = ring.GetPointCount()
            for p in range(numpoints):
                    lon, lat, z = ring.GetPoint(p)
                    pointsX.append(lon)
                    pointsY.append(lat)
            count += 1
    elif (geom.GetGeometryName() == 'POLYGON'):
        ring = geom.GetGeometryRef(0)
        numpoints = ring.GetPointCount()
        pointsX = []; pointsY = []
        for p in range(numpoints):
                lon, lat, z = ring.GetPoint(p)
                pointsX.append(lon)
                pointsY.append(lat)    
    else:
        sys.exit("ERROR: Geometry needs to be either Polygon or Multipolygon")    
    xmin = min(pointsX)
    xmax = max(pointsX)
    ymin = min(pointsY)
    ymax = max(pointsY)    
    # Specify offset and rows and columns to read
    xoff = int((xmin - xOrigin)/pixelWidth)
    #if ind > 6:
    #    xoff = xoff + 1
    yoff = int((yOrigin - ymax)/pixelWidth)
    #if yoff > 0:
    #    yoff = yoff + 1
    xcount = int((xmax - xmin)/pixelWidth)+1
    ycount = int((ymax - ymin)/pixelWidth)+1
#    with open("/scratch/pj276/ifl_corridors/geodata/hires_global_v1_0_2013/diagnostic_log.txt", "a") as text_file:
#        text_file.write("xoff = " + str(xoff) + "\n")
#        text_file.write("yoff = " + str(yoff) + "\n")
#        text_file.write("xorigin = " + str(xOrigin) + "\n")
#        text_file.write("yorigin = " + str(yOrigin) + "\n")
#        text_file.write("xcount = " + str(xcount) + "\n")
#        text_file.write("ycount = " + str(ycount) + "\n")
    # Create memory target raster
    driver = gdal.GetDriverByName('GTiff')
    target_ds = driver.Create(newRasterfn, xcount, ycount, 1, gdal.GDT_Byte)
    target_ds.SetGeoTransform((
            xmin, pixelWidth, 0,
            ymax, 0, pixelHeight,
    ))
    # Create for target raster the same projection as for the value raster
    raster_srs = osr.SpatialReference()
    raster_srs.ImportFromWkt(raster.GetProjectionRef())
    target_ds.SetProjection(raster_srs.ExportToWkt())
    # Read raster as array
    dataraster = raster.ReadAsArray(xoff, yoff, xcount, ycount).astype(np.int8)
    # Take max or sum of dataraster across vertical dimension if it has one
    # Then write array to raster.
    #mdata = np.ma.array(dataraster, mask=(dataraster == 0))
    if len(dataraster.shape) == 3:
        if mosfun == 'max':
            drmax = dataraster.max(axis=0)
        elif mosfun == 'sum':
            drmax = dataraster.sum(axis=0)
        else:
            drmax = dataraster.max(axis=0)
        tdsb = target_ds.GetRasterBand(1)
        tdsb.WriteArray(drmax)
        tdsb.FlushCache()
        drmax = None
        tdsb = None
        target_ds = None
        raster = None
        dataraster = None
        shp = None
    else:
        tdsb = target_ds.GetRasterBand(1)
        tdsb.WriteArray(dataraster)
        tdsb.FlushCache()
        tdsb = None
        target_ds = None
        raster = None
        dataraster = None
        shp = None

# This function takes a multiband input vrt and makes a mosaic
# of the vrt. Takes either 'max' or 'sum' as a mosaic function (mosfun).
# gdtype is gdal pixel type object e.g. gdal.GDT_Float32
# nptype is corresponding numpy pixel type as string e.g. 'float32'
def tile_shingle(input_value_raster, newRasterfn, gdtype, nptype, mosfun='max'):
    # Open data
    raster = gdal.Open(input_value_raster)
    if raster is not None:
        # Get raster georeference info
        transform = raster.GetGeoTransform()
        # Create memory target raster (same as input)
        driver = gdal.GetDriverByName('GTiff')
        target_ds = driver.Create(newRasterfn, raster.RasterXSize, raster.RasterYSize, 1, gdtype, ['COMPRESS=LZW'])
        target_ds.SetGeoTransform(transform)
        # Create for target raster the same projection as for the value raster
        raster_srs = osr.SpatialReference()
        raster_srs.ImportFromWkt(raster.GetProjectionRef())
        target_ds.SetProjection(raster_srs.ExportToWkt())
        # Read raster to array
        dataraster = raster.ReadAsArray().astype(nptype)
        # Convert nan to 0, infinite to large numbers
        dataraster = np.nan_to_num(dataraster)
        # Take max of shingled vrt across vertical dimension if it has one
        # Then write array to raster.
        #mdata = np.ma.array(dataraster, mask=(dataraster == 0))
        if len(dataraster.shape) == 3:
            if mosfun == 'sum':
                drmax = dataraster.sum(axis=0)
            elif mosfun == 'mean':
                drmax = dataraster[dataraster > 0].mean(axis=0)
            else:
                drmax = dataraster.max(axis=0)
            if drmax.max() > 0:
                tdsb = target_ds.GetRasterBand(1)
                tdsb.WriteArray(drmax)
                tdsb.FlushCache()
                drmax = None
                tdsb = None
                target_ds = None
                raster = None
                dataraster = None
            else:
                drmax = None
                target_ds = None
                raster = None
                dataraster = None
                os.remove(newRasterfn)
        else:
            tdsb = target_ds.GetRasterBand(1)
            tdsb.WriteArray(dataraster)
            tdsb.FlushCache()
            tdsb = None
            target_ds = None
            raster = None
            dataraster = None

# Alternative randomized shortest path function. Seems to work a bit
# faster than the first one I wrote. rfn is a sparse transition matrix.
# index 1 and index 2 are indices (cell numbers I think) that give the
# start and end points of the corridor. Theta is level of path randomness.
# theta = 0.5 is a good starting point
def rsc2(rfn,index1,index2,theta,nrow,ncol):
    trR = rfn
    trR.data = 1/trR.data
    nr = trR.shape[0] # Get number of transition matrix rows
    Id = identity(nr, dtype='float32', format='dia')# Create sparse identity matrix
    rs = csr_matrix(trR.sum(1)) # Get row sums of transition matrix as sparse matrix
    rs[rs>0] = 1/rs[rs>0] # Take inverse of row sums greater than 0. Should still be sparse
    #rs.data = 1/rs.data
    P = trR.multiply(rs) # Maybe this is a transition probability matrix. Multiply outputs csr so convert back to csc
    W = trR
    W.data = np.exp(-theta * trR.data)
    W = W.multiply(P)
    #nc = trR.shape[0] # Get number of cells in raster used to derive transition matrix (shortcut is to grab one of the dimensions of the transition matrix)
    Ij = identity(nr, dtype='float32', format='dia') # Identity matrix
    # Indices. Be careful of zero indexing
    ci = np.matrix(index1)
    cj = np.matrix(index2)
    Ij.data[:,cj] = 1 - 1 / float(len(cj))
    Wj = Ij * W
    ei = [0]*nr
    ei[int(ci)] = 1 / float(len(ci))
    ej = [0]*nr
    ej[int(cj)] = 1 / float(len(cj))
    IdMinusWj = Id - Wj
    zci = spl.spsolve(IdMinusWj.transpose(),csr_matrix(ei,dtype="float32").transpose(),use_umfpack=True)
    zcj = spl.spsolve(IdMinusWj,csr_matrix(ej,dtype="float32").transpose(),use_umfpack=True)
    zcij = sum(np.array(ei)*zcj)
    N = diags(zci) * Wj * diags(zcj) / zcij
    # TOTAL number of passages?
    # Get maximum of row sums and column sums
    # n = np.maximum(np.sum(N,1,dtype='float32'), np.sum(N,0,dtype='float32').T)
    # n = np.squeeze(np.array(n)) # Get rid of extra dimension
    # NET number of passages?
    # Get skew symmetric matrix part
    nNet = abs((N-N.T)/np.float(2))
    nn = np.maximum(nNet.sum(1), nNet.sum(0).T)
    nn = np.squeeze(np.array(nn)) # Get rid of extra dimension
    # Set origin and destination cell values
    nn[[int(ci),int(cj)]] = 2 * nn[[int(ci),int(cj)]]
    nn = nn.reshape(nrow,ncol)
    return(nn)

# This function works but is not set up to run without moving arguments
# to the function call. It does the same thing (I think) as rsc2
def rsc(rfn,index1,index2,theta,nrow,ncol):
    #theta = 0.5
    #trR = triu(csc_matrix(trb, dtype=np.float32),format="csc") # return upper triangle of transition matrix as sparse, compressed, column oriented format
    trR = rfn
    trR.data = 1/trR.data
    nr = trR.shape[0] # Get number of transition matrix rows
    Id = identity(nr, dtype='float32', format='dia')# Create sparse identity matrix
    rs = csr_matrix(csr_matrix.sum(trR,1),dtype=np.float32) # Get row sums of transition matrix as sparse matrix
    rs[rs>0] = 1/rs[rs>0] # Take inverse of row sums greater than 0. Should still be sparse
    P = csr_matrix(trR.multiply(rs)) # Maybe this is a transition probability matrix. Multiply outputs csr so convert back to csc
    W = trR
    aa = trR.data
    W.data = np.exp(-theta * aa)
    W = W.multiply(P)
    #nc = trR.shape[0] # Get number of cells in raster used to derive transition matrix (shortcut is to grab one of the dimensions of the transition matrix)
    Ij = identity(nr, dtype='float32', format='dia') # Identity matrix
    # Hard code indices for testing. Be careful of zero indexing
    #ci = np.matrix(2)
    #cj = np.matrix(8)
    ci = np.matrix(index1)
    cj = np.matrix(index2)
    Ij.data[:,cj] = 1 - 1 / len(cj)
    #Wj = Ij.dot(W) # Matrix multiply
    Wj = Ij * W
    Wj = csr_matrix(Wj)
    ei = [0]*nr
    ei[int(ci)] = 1 / len(ci)
    ej = [0]*nr
    ej[int(cj)] = 1 / len(cj)
    IdMinusWj = Id - Wj
    IdMinusWj = csr_matrix(IdMinusWj)
    zci = spl.spsolve(np.transpose(IdMinusWj),csr_matrix(np.transpose(np.matrix(ei)),dtype="float32"))
    zcj = spl.spsolve(IdMinusWj,csr_matrix(np.transpose(np.matrix(ej)),dtype="float32"))
    zcij = sum(np.array(ei)*zcj)
    N = csr_matrix(np.diag(zci)) * Wj * csr_matrix(np.diag(zcj)) / zcij
    # Convert to dense matrix format
    N = N.todense()
    # TOTAL number of passages?
    # Get maximum of row sums and column sums
    # n = np.maximum(np.sum(N,1,dtype='float32'), np.sum(N,0,dtype='float32').T)
    # n = np.squeeze(np.array(n)) # Get rid of extra dimension
    # NET number of passages?
    # Get skew symmetric matrix part
    nNet = np.abs((N-N.T)/2)
    # Get maximum of row sums and column sums
    nn = np.maximum(np.sum(nNet,1,dtype='float32'), np.sum(nNet,0,dtype='float32').T)
    nn = np.squeeze(np.array(nn)) # Get rid of extra dimension
    # Set origin and destination cell values
    nn[[int(ci),int(cj)]] = 2 * nn[[int(ci),int(cj)]]
    nn = nn.reshape(nrow,ncol)
    return(nn)

# Minor modification from 
# https://stackoverflow.com/questions/30199070/how-to-create-a-4-or-8-connected-adjacency-matrix
# Diagonal weights set to sqrt(2). Orthogonal weights set to 1.
# Calculation for conductance is then conductance value/cellres*weight
def connected_adjacency(image, connect, patch_size=(1, 1)):
    """
    Creates an adjacency matrix from an image where nodes are considered adjacent 
    based on 4-connected or 8-connected pixel neighborhoods.
    :param image: 2 or 3 dim array
    :param connect: string, either '4' or '8'
    :param patch_size: tuple (n,m) used if the image will be decomposed into 
                   contiguous, non-overlapping patches of size n x m. The 
                   adjacency matrix will be formed from the smaller sized array
                   e.g. original image size = 256 x 256, patch_size=(8, 8), 
                   then the image under consideration is of size 32 x 32 and 
                   the adjacency matrix will be of size 
                   32**2 x 32**2 = 1024 x 1024
    :return: adjacency matrix as a sparse matrix (type=scipy.sparse.csr.csr_matrix)
    """
    r, c = image.shape[:2]
    r = r / patch_size[0]
    c = c / patch_size[1]
    if connect == '4':
        # constructed from 2 diagonals above the main diagonal
        d1 = np.tile(np.append(np.ones(c-1), [0]), r)[:-1]
        d2 = np.ones(c*(r-1))
        upper_diags = diags([d1, d2], [1, c])
        return upper_diags + upper_diags.T
    elif connect == '8':
        # constructed from 4 diagonals above the main diagonal
        d1 = np.tile(np.append(np.ones(c-1), [0]), r)[:-1]
        d2 = np.append([0], d1[:c*(r-1)])
        d3 = np.ones(c*(r-1))
        d4 = d2[1:-1]
        d4[d4==1] = 2.0**0.5
        upper_diags = diags([d1, d2, d3, d4], [1, c-1, c, c+1])
        return upper_diags + upper_diags.T
    else:
        raise ValueError('Invalid parameter \'connect\'={connect}, must be "4" or "8".'
                     .format(connect=repr(connect)))

# Function to get column number from cell number (based on 1 indexing but converted to python indexing)
# matcols is the number of columns in the array
def cell2inds(cellnum, matcols):
    cellnum = cellnum+1
    acol = cellnum-(math.ceil(cellnum/float(matcols))*float(matcols)-(float(matcols)-1))+1
    arow = math.ceil(cellnum/float(matcols))
    acol = int(acol-1)
    arow = int(arow-1)
    return((arow,acol))

# Return cell number from row,column indices. I believe cell numbering
# runs from upper left to right and down.
# matcols is the number of columns in the array.
# give row and column indices using python indexing
# Does not error check for compatible dimensions yet
def inds2cell(arow,acol,matcols):
    arow=arow+1
    acol=acol+1
    xx = ((arow-1)*matcols)+1+acol-1
    return(xx-1) # subtract 1 for python indexing

# This function uses morphological image processing to identify
# patch pixels given a structuring element and a binary raster.
# Input array should have target class coded as 1, all else zero.
# givecore is 0 or 1 and specifies whether the core pixel array
# should be returned along with the patch pixel map.
def idpatchpix(inarray, st_elem4, st_elem8, givecore=1):
    # Erosion
    ero1 = binary_erosion(inarray, st_elem8)
    # Dilations
    # First dilation
    dil1 = binary_dilation(ero1,st_elem4)
    # Add orthogonally connected pixels to eroded image
    newpix = np.where((dil1==1) & (inarray==1), 1, 0)
    # Get pixel baseline count
    pc = np.count_nonzero(newpix)
    # Second dilation
    dil1 = binary_dilation(newpix,st_elem4)
    # Add orthogonal pixels
    newpix = np.where((dil1==1) & (inarray==1), 1, 0)
    # Get difference between first and second orthogonal pixel addition
    cdiff = np.count_nonzero(newpix)-pc
    # Keep dilating until no more orthogonally connected pixels are added
    while cdiff > 0:
        pc = np.count_nonzero(newpix)    
        dil1 = binary_dilation(newpix,st_elem4)
        newpix = np.where((dil1==1) & (inarray==1), 1, 0)
        cdiff = np.count_nonzero(newpix)-pc
    # Identify forest patches
    patchpix = np.where((inarray == 1) & (newpix == 0), 1, 0)
    # Identify other forest
    of = np.where((patchpix == 0) & (inarray ==1 ), 1, 0)
    # Return only patches if givecore = 0.
    if givecore == 0:
        return(patchpix)
    else:
        return(patchpix, ero1, of)

# This function uses morphological image processing to identify
# patch pixels given a structuring element and a binary raster.
# Input array should have target class coded as 1, all else zero.
# givecore is 0 or 1 and specifies whether the core pixel array
# should be returned along with the patch pixel map.
def idpatchpixbool(inarray, st_elem4, st_elem8, givecore=1):
    # Erosion
    ero1 = mahotas.morph.erode(inarray, st_elem8)
    # Dilations
    # First dilation
    dil1 = mahotas.morph.dilate(ero1,st_elem4)
    # Add orthogonally connected pixels to eroded image
    newpix = np.where((dil1==True) & (inarray==True), True, False)
    # Get pixel baseline count
    pc = np.count_nonzero(newpix)
    # Second dilation
    dil1 = mahotas.morph.dilate(newpix,st_elem4)
    # Add orthogonal pixels
    newpix = np.where((dil1==True) & (inarray==True), True, False)
    # Get difference between first and second orthogonal pixel addition
    cdiff = np.count_nonzero(newpix)-pc
    # Keep dilating until no more orthogonally connected pixels are added
    while cdiff > 0:
        pc = np.count_nonzero(newpix)    
        dil1 = mahotas.morph.dilate(newpix,st_elem4)
        newpix = np.where((dil1==True) & (inarray==True), True, False)
        cdiff = np.count_nonzero(newpix)-pc
    # Identify forest patches
    patchpix = np.where((inarray == True) & (newpix == False), True, False)
    # Identify other forest
    of = np.where((patchpix == False) & (inarray == True), True, False)
    # Return only patches if givecore = 0.
    if givecore == 0:
        return(patchpix)
    else:
        return(patchpix, ero1, of)

# This function fills in pixels that are both core edge and perforation
# edge. It gives core edge priority, reclassifying perforation edge as
# core edge.
def buildnewedge(inperf, inedge):
    # Dilations
    # First dilation (use default structuring element, 3x3 cross)
    dil1 = binary_dilation(inedge)
    # Add orthogonally connected pixels to eroded image
    newpix = np.where((dil1==1) & (inperf==1), 1, 0)
    newpix = newpix.astype('int8')
    # Get pixel baseline count
    pc = np.count_nonzero(newpix)
    # Second dilation
    dil1 = binary_dilation(newpix)
    # Add orthogonal pixels
    newpix = np.where((dil1==1) & (inperf==1), 1, 0)
    newpix = newpix.astype('int8')   
    # Get difference between first and second orthogonal pixel addition
    cdiff = np.count_nonzero(newpix)-pc
    # Keep dilating until no more orthogonally connected pixels are added
    while cdiff > 0:
        pc = np.count_nonzero(newpix)    
        dil1 = binary_dilation(newpix)
        newpix = np.where((dil1==1) & (inperf==1), 1, 0)
        newpix = newpix.astype('int8')
        cdiff = np.count_nonzero(newpix)-pc
    return(newpix)

# This function fills in pixels that are both core edge and perforation
# edge. It gives core edge priority, reclassifying perforation edge as
# core edge.
def buildnewedgebool(inperf, inedge):
    # Dilations
    # First dilation (use default structuring element, 3x3 cross)
    dil1 = mahotas.morph.dilate(inedge)
    # Add orthogonally connected pixels to eroded image
    newpix = np.where((dil1==True) & (inperf==True), True, False)
    #newpix = newpix.astype('int8')
    # Get pixel baseline count
    pc = np.count_nonzero(newpix)
    # Second dilation
    dil1 = mahotas.morph.dilate(newpix)
    # Add orthogonal pixels
    newpix = np.where((dil1==True) & (inperf==True), True, False)
    #newpix = newpix.astype('int8')   
    # Get difference between first and second orthogonal pixel addition
    cdiff = np.count_nonzero(newpix)-pc
    # Keep dilating until no more orthogonally connected pixels are added
    while cdiff > 0:
        pc = np.count_nonzero(newpix)    
        dil1 = mahotas.morph.dilate(newpix)
        newpix = np.where((dil1==True) & (inperf==True), True, False)
        #newpix = newpix.astype('int8')
        cdiff = np.count_nonzero(newpix)-pc
    return(newpix)

# This function takes an input array and a structuring element and clips
# the array using the dimensions of the structuring element. Used
# to get rid of padding from morphological operations.
def trimarray(inarray, selem):
    inarray = inarray[0+selem*2:inarray.shape[0]-selem*2,0+selem*2:inarray.shape[1]-selem*2]
    return(inarray)
# Reclassify forest integrity to guidos compliant
# inval is the threshold value below which the raster is coded to background
# ndval is the value in the input raster that is NoData
# rclzero is the value to reclass zero to. Setting it at 0 means no reclass of 
# zero will take place.
def rclguidos(inarray, inval, ndval, rclzero):
    if rclzero > 0:
        inarray = np.where(inarray == 0, rclzero, inarray)
    inarray = np.where(inarray == ndval, 0, inarray)
    inarray = np.where((inarray < inval) & (inarray != 0), 1, inarray)
    inarray = np.where((inarray >= inval) & (inarray != 0) & (inarray != 1), 2, inarray)
    inarray = inarray.astype('uint8')
    return(inarray)

# Function takes a shapefile name (full path, with extension), an index, and x and y intervals as integers.
# Also takes an output shapefile path (with extension).
# Returns the name of the output point shapefile
# Creates a temporary polygon shapefile in the directory of the input shapefile
# inval should be an integer
def regular_points(inshp, infield, inval, xint, yint, dosnap=0, snap_x=1, snap_y=1):
    # Range function for floats
    def floatrange(start, stop, step):
        while start < stop:
            yield start
            start += step
    # Create a reader instance
    r = shp.Reader(inshp)
    fnames = [i[0] for i in r.fields]
    findex = fnames.index(infield)-1 # need to subtract 1 to account for deletion field
    # Create a writer instance
    w = shp.Writer(shapeType=shp.POLYGON)
    # Copy the fields to the writer
    w.fields = list(r.fields)
    # Grab the geometry and records from all features 
    # with the correct gridcode 
    selection = [] 
    for rec in enumerate(r.records()):
       if rec[1][findex] == inval:
          selection.append(rec) 
    # Add the geometry and records to the writer
    for rec in selection:
       w._shapes.append(r.shape(rec[0]))
       w.records.append(rec[1])
    # Save the new shapefile
    tpoly = inshp.replace('.shp','_' + str(rec[0]) + '_' + 'poly.shp')
    w.save(tpoly) 
    ashape = r.shape(rec[0])
    xpts = [] # list to hold x coords
    ypts = [] # list to hold y coords
    # Populate lists of x and y coords
    for u in ashape.points:
        xpts.append(u[0])
        ypts.append(u[1])
    # Get mins and maxes
    xmin = min(xpts)
    xmax = max(xpts)
    ymin = min(ypts)
    ymax = max(ypts)
    # Snap if requested
    if dosnap == 1:
        # Get number of columns and rows from extent, plus 1
        newxres = int((xmax - xmin) / xint) + 1
        newyres = int((ymax - ymin) / yint) + 1
        # Adjust xmin and ymax to align with the snap grid
        #newymax = (snap_y-(math.ceil(abs((ymax-snap_y)/30))*30))+30*5
        #newxmin = ((math.ceil(abs((xmin-snap_x)/30))*30)+snap_x)-30*5
        if ymax < snap_y:
            newymax = (snap_y-(math.ceil(abs((ymax-snap_y)/yint))*yint))+yint*1
        else:
            newymax = (snap_y+(math.ceil(abs((ymax-snap_y)/yint))*yint))+yint*1
        if xmin > snap_x:
            newxmin = (snap_x+(math.ceil(abs((xmin-snap_x)/xint))*xint))-xint*1
        else:
            newxmin = (snap_x-(math.ceil(abs((xmin-snap_x)/xint))*xint))-xint*1  
        newymin = newymax - newyres*yint
        newxmax = newxmin + newxres*xint
        ll = [newxmin,newymin]
        ur = [newxmax,newymax]
    else:
        ll = [xmin,ymin] # lower left coordinate pair of polygon's extent
        ur = [xmax,ymax] # upper right
    # Get polygon as wkt
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(tpoly, 0)
    layer = dataSource.GetLayer()
    for feature in layer:
        geom = feature.GetGeometryRef()
    x_interval = xint # define interval
    y_interval = yint # define interval
    # Output coordinate list
    xoco = []
    yoco = []
    # Loop through and create points
    for xco in floatrange(ll[0], ur[0], x_interval):
        for yco in floatrange(ll[1], ur[1], y_interval):
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(xco,yco)
            if point.Within(geom):
                xoco.append(xco)
                yoco.append(yco)    
    point = None
    plocations = zip(xoco,yoco)
    # Clean up
    layer = None
    geom = None
    if os.path.exists(tpoly):
        shpDriver = ogr.GetDriverByName("ESRI Shapefile")
        shpDriver.DeleteDataSource(tpoly)
        shpDriver = None
    return(plocations)

# Function takes a shapefile name (full path, with extension), an index, and x and y intervals as integers.
# Also takes an output shapefile path (with extension).
# Also takes x,y coordinates of a centroid or representative point of a shapefile
# Creates a temporary polygon shapefile in the directory of the input shapefile
# inval should be an integer
#('C:\Users\pj276\Google Drive\Projects\ifl_corridors\outputs\hland_tropics_ssu_sa_subset2\hland_tropics_ssu_sa_subset2.shp','GRIDCODE',2539,10000,10000,1356579.82, -1359995.12)
def regular_points_from_center(inshp, infield, inval, xint, yint):
    # Range function for floats
    def floatrange(start, stop, step):
        while start < stop:
            yield start
            start += step
    # If final point shapefile doesn't exist already, create it.
    # Create a reader instance
    r = shp.Reader(inshp)
    fnames = [i[0] for i in r.fields]
    findex = fnames.index(infield)-1 # need to subtract 1 to account for deletion field
    # Create a writer instance
    w = shp.Writer(shapeType=shp.POLYGON)
    # Copy the fields to the writer
    w.fields = list(r.fields)
    # Grab the geometry and records from all features 
    # with the correct gridcode 
    selection = [] 
    for rec in enumerate(r.records()):
       if rec[1][findex] == inval:
          selection.append(rec) 
    # Add the geometry and records to the writer
    for rec in selection:
       w._shapes.append(r.shape(rec[0]))
       w.records.append(rec[1])
    # Save the new shapefile
    tpoly = inshp.replace('.shp','_' + str(rec[0]) + '_' + 'poly.shp')
    w.save(tpoly) 
    ashape = r.shape(rec[0])
    xpts = [] # list to hold x coords
    ypts = [] # list to hold y coords
    # Populate lists of x and y coords
    for u in ashape.points:
        xpts.append(u[0])
        ypts.append(u[1])
    # Use fiona and geopandas to dissolve features if there are multiple ones with the same gridcode
    reader = fiona.open(tpoly)
    myproj = reader.crs
    xx = gpd.GeoDataFrame.from_features((x for x in reader if x['properties'][infield]==inval))
    xx = xx.dissolve(by=infield, aggfunc='sum')
    xx.crs = myproj
    xx[infield] = inval
    # Get mins and maxes
    xmin = min(xpts)
    xmax = max(xpts)
    ymin = min(ypts)
    ymax = max(ypts)
    # Get representative point (need to use the inval as a dict key)
    rp = xx.representative_point()
    cx = rp.x.get(inval)
    cy = rp.y.get(inval)
    # Get new lower left and upper right based on spacing from center point
    newll = [cx - (math.ceil((cx-xmin)/xint)*xint), cy - (math.ceil((cy-ymin)/yint)*yint)]
    newur = [cx + (math.ceil((xmax-cx)/xint)*xint), cy + (math.ceil((ymax-cy)/yint)*yint)]
    # Get polygon as wkt
    driver = ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.Open(tpoly, 0)
    layer = dataSource.GetLayer()
    for feature in layer:
        geom = feature.GetGeometryRef()
    # Output coordinate list
    xoco = []
    yoco = []
    # Loop through and create points
    for xco in floatrange(newll[0], newur[0], xint):
        for yco in floatrange(newll[1], newur[1], yint):
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(xco,yco)
            xoco.append(xco)
            yoco.append(yco)
##            if point.Within(geom):
##                xoco.append(xco)
##                yoco.append(yco)    
    point = None
    plocations = zip(xoco,yoco)
    # Clean up
    layer = None
    geom = None
    if os.path.exists(tpoly):
        shpDriver = ogr.GetDriverByName("ESRI Shapefile")
        shpDriver.DeleteDataSource(tpoly)
        shpDriver = None
    return(plocations)

# Return x,y of the representative point of a shapefile
def representative_point(inshp, infield, inval):
    # Create a reader instance
    r = shp.Reader(inshp)
    fnames = [i[0] for i in r.fields]
    findex = fnames.index(infield)-1 # need to subtract 1 to account for deletion field
    # Create a writer instance
    w = shp.Writer(shapeType=shp.POLYGON)
    # Copy the fields to the writer
    w.fields = list(r.fields)
    # Grab the geometry and records from all features 
    # with the correct gridcode 
    selection = [] 
    for rec in enumerate(r.records()):
       if rec[1][findex] == inval:
          selection.append(rec) 
    # Add the geometry and records to the writer
    for rec in selection:
       w._shapes.append(r.shape(rec[0]))
       w.records.append(rec[1])
    # Save the new shapefile
    tpoly = inshp.replace('.shp','_' + str(rec[0]) + '_' + 'poly.shp')
    w.save(tpoly) 
    # Use fiona and geopandas to dissolve features if there are multiple ones with the same gridcode
    reader = fiona.open(tpoly)
    myproj = reader.crs
    xx = gpd.GeoDataFrame.from_features((x for x in reader if x['properties'][infield]==inval))
    xx = xx.dissolve(by=infield, aggfunc='sum')
    xx.crs = myproj
    xx[infield] = inval
    # Get representative point (need to use the inval as a dict key)
    rp = xx.representative_point()
    cx = rp.x.get(inval)
    cy = rp.y.get(inval)
    plocations = [(cx,cy)] # This needs to be a list of a tuple for the coord to point function
    # Clean up
    if os.path.exists(tpoly):
        shpDriver = ogr.GetDriverByName("ESRI Shapefile")
        shpDriver.DeleteDataSource(tpoly)
        shpDriver = None
    return(plocations)

# This function finds the envelope of the specified polygon in a shapefile
# Specification is based on a field value. Returns a list of xmin, xmax, ymin, and ymax.
# inval should be an integer
def poly_envelope(inshp, infield, inval):
    # Range function for floats
    def floatrange(start, stop, step):
        while start < stop:
            yield start
            start += step
    # Create a reader instance
    r = shp.Reader(inshp)
    fnames = [i[0] for i in r.fields]
    findex = fnames.index(infield)-1 # need to subtract 1 to account for deletion field
    # Grab the geometry and records from all features 
    # with the correct gridcode 
    selection = [] 
    for rec in enumerate(r.records()):
       if rec[1][findex] == inval:
          selection.append(rec)
    bblist = [] # get list of bounding boxes (lower left x,y and upper right x,y)
    for rec in selection:
        ashape = r.shape(rec[0])
        bblist.append(ashape.bbox)
    xmins = [i[0] for i in bblist]
    ymins = [i[1] for i in bblist]
    xmaxs = [i[2] for i in bblist]
    ymaxs = [i[3] for i in bblist]
    xmin = min(xmins)
    xmax = max(xmaxs)
    ymin = min(ymins)
    ymax = max(ymaxs)
    return([xmin,xmax,ymin,ymax])
                   
# Merge shapefiles
# Merge a bunch of shapefiles with attributes quickly!
# Cribbed from here http://geospatialpython.com/2011/02/merging-lots-of-shapefiles-quickly.html
def point_merge(inlist, outname):
    w = shp.Writer()
    for f in inlist:
        r = shp.Reader(f)
        w._shapes.extend(r.shapes())
        w.records.extend(r.records())
    w.fields = list(r.fields)
    w.save(outname)

# Add numeric field to shapefile
# Integer field type is 'N'
# Can check here for more info https://pypi.python.org/pypi/pyshp
# Cribbed from here http://geospatialpython.com/2013/04/add-field-to-existing-shapefile.html
# Example - add_shp_field(opath, 'GCODE', 'N', 8, 431, opath)
# I think if outfile is same is infile, it overwrites original
def add_shp_field(inshape, infield, fieldtype, fieldwidth, fillval, outshape):
    # Read in our existing shapefile
    r = shp.Reader(inshape)
    # Create a new shapefile in memory
    w = shp.Writer()
    # Copy over the existing fields
    w.fields = list(r.fields)
    # Add our new field using the pyshp API
    w.field(infield, fieldtype, fieldwidth)
    # Loop through each record, add a column and data.
    # you could also just
    # insert a blank string or NULL DATA number
    # as a place holder
    for rec in r.records():
        rec.append(fillval)
        # Add the modified record to the new shapefile 
        w.records.append(rec)
    # Copy over the geometry without any changes
    w._shapes.extend(r.shapes())
    # If creating a new shapefile, copy over projection file
    if outshape != inshape:
        # Read in template proj
        prjref = inshape.replace(".shp",".prj")
        f = open(prjref)
        fr = f.read()
        f = None
        # Replace shp extension with prj extension for output file
        prjname = outshape.replace(".shp",".prj")
        prj = open(prjname, "w")
        prj.write(fr)
        prj.close()
    # Save as a new shapefile (or write over the old one)
    w.save(outshape)

# This is the rSPDist function from gdistance. Need to convert to python.
#function (tr, ci, cj, theta, totalNet, method) 
#{
#    trR <- tr
#    trR@x <- 1/trR@x
#    nr <- dim(tr)[1]
#    Id <- Diagonal(nr)
#    P <- .normalize(tr, "row")
#    if (method == 1) {
#        W <- trR
#        W@x <- exp(-theta * trR@x)
#        W <- W * P
#    }
#    else {
#        adj <- adjacencyFromTransition(tr)
#        W <- trR
#        W[adj] <- exp(-theta * -log(P[adj]))
#    }
#    D <- matrix(0, nrow = length(ci), ncol = length(cj))
#    for (j in 1:length(cj)) {
#        Ij <- Diagonal(nr)
#        Ij[cj[j], cj[j]] <- 0
#        Wj <- Ij %*% W
#        IdMinusWj <- as((Id - Wj), "dgCMatrix")
#        ej <- rep(0, times = nr)
#        ej[cj[j]] <- 1
#        zcj <- solve(IdMinusWj, ej)
#        for (i in 1:length(ci)) {
#            ei <- rep(0, times = nr)
#            ei[ci[i]] <- 1
#            zci <- solve(t(IdMinusWj), ei)
#            zcij <- sum(ei * zcj)
#            N <- (Diagonal(nr, as.vector(zci)) %*% Wj %*% Diagonal(nr, 
#                as.vector(zcj)))/zcij
#            if (totalNet == "net") {
#                N <- skewpart(N) * 2
#                N@x[N@x < 0] <- 0
#            }
#            D[i, j] <- sum(trR * N)
#        }
#    }
#    return(D)
#}

# This function creates a convex hull from a shapefile. It gives the feature
# a field called 'id' and gives it a value of 1. inShapefile needs a .shp extension.
# So does outShapefile. Give buffer in meters if you want to buffer the hull.
# Adapted from https://pcjericks.github.io/py-gdalogr-cookbook/vector_layers.html#save-the-convex-hull-of-all-geometry-from-an-input-layer-to-an-output-layer
def convexhull(inShapefile, outShapefile, buff=None):
    # Get a Layer
    inDriver = ogr.GetDriverByName("ESRI Shapefile")
    inDataSource = inDriver.Open(inShapefile, 0)
    inLayer = inDataSource.GetLayer()
    # Collect all Geometry
    geomcol = ogr.Geometry(ogr.wkbGeometryCollection)
    for feature in inLayer:
        geomcol.AddGeometry(feature.GetGeometryRef())
    # Calculate convex hull
    convexhull = geomcol.ConvexHull()
    # Buffer if desired
    if buff is not None:
        convexhull = convexhull.Buffer(buff)
    # Save extent to a new Shapefile
    outDriver = ogr.GetDriverByName("ESRI Shapefile")
    # Remove output shapefile if it already exists
    if os.path.exists(outShapefile):
        outDriver.DeleteDataSource(outShapefile)
    # Create the output shapefile
    outDataSource = outDriver.CreateDataSource(outShapefile)
    outLayer = outDataSource.CreateLayer(os.path.splitext(outShapefile)[0], geom_type=ogr.wkbPolygon)
    # Add an ID field
    idField = ogr.FieldDefn("id", ogr.OFTInteger)
    outLayer.CreateField(idField)
    # Create the feature and set values
    featureDefn = outLayer.GetLayerDefn()
    feature = ogr.Feature(featureDefn)
    feature.SetGeometry(convexhull)
    feature.SetField("id", 1)
    outLayer.CreateFeature(feature)
    feature = None
    # Read in template proj
    prjref = inShapefile.replace(".shp",".prj")
    f = open(prjref)
    fr = f.read()
    f = None
    # Replace shp extension with prj extension for output file
    prjname = outShapefile.replace(".shp",".prj")
    prj = open(prjname, "w")
    prj.write(fr)
    prj.close()
    # Save and close DataSource
    inDataSource = None
    outDataSource = None

def shp2srs(inshp):
    prjf = inshp.replace('.shp','.prj') # Get projection file name
    f = open(prjf) # Open projection file
    fr = f.read() # Read to variable
    srs = osr.SpatialReference(wkt=fr) # Convert to osr spatial reference
    f = None # Close file
    return srs

def vrtflist(invrt):
    raster = gdal.Open(invrt)
    flist = raster.GetFileList()
    return flist

# Adapted from
# https://gis.stackexchange.com/questions/178765/intersecting-two-shapefiles-from-python-or-command-line
# I think the arguments should be strings. Output name should have .shp extension.
# This returns a list of dictionaries where the keys are 'field1' and 'field2'.
# Iterate through and use dict.get() to return the field of interest.
# The function expects a 'uid' field. It will add it if it doesn't exist
def polinpol(shp1, shp2, field1='uid', field2='uid'):
    g1 = gpd.GeoDataFrame.from_file(shp1)
    if 'uid' not in list(g1):
        g1['uid'] = None
        for index, row in g1.iterrows():
            g1.loc[index, 'uid'] = index
    g2 = gpd.GeoDataFrame.from_file(shp2)
    if 'uid' not in list(g2):
        g2['uid'] = None
        for index, row in g2.iterrows():
            g2.loc[index, 'uid'] = index
    data = []
    for index, s1 in g1.iterrows():
        for index2, s2 in g2.iterrows():
           if s1['geometry'].intersects(s2['geometry']):
              #data.append({'geometry': s1['geometry'].intersection(s2['geometry']), field1:s1[field1], field2: s2[field2], 'area':s1['geometry'].intersection(s2['geometry']).area})
              #data.append({field1:s1[field1], field2: s2[field2]})
              data.append({'field1':s1[field1], 'field2': s2[field2]})
    #df = gpd.GeoDataFrame(data,columns=['geometry', field1, field2,'area'])
    #return df
    return data
    #df.to_file(outsfname)

# Write intersection shapefile to disk
# The function expects a 'uid' field. It will add it if it doesn't exist
def intersectionGeom(shp1, shp2, field1='uid', field2='uid', write='no', sfname='default'):
    g1 = gpd.GeoDataFrame.from_file(shp1)
    if 'uid' not in list(g1):
        g1['uid'] = None
        for index, row in g1.iterrows():
            g1.loc[index, 'uid'] = index
    g2 = gpd.GeoDataFrame.from_file(shp2)
    if 'uid' not in list(g2):
        g2['uid'] = None
        for index, row in g2.iterrows():
            g2.loc[index, 'uid'] = index
    data = []
    for index, s1 in g1.iterrows():
        for index2, s2 in g2.iterrows():
           if s1['geometry'].intersects(s2['geometry']):
              data.append({'geometry': s1['geometry'].intersection(s2['geometry']), field1:s1[field1], 'area':s1['geometry'].intersection(s2['geometry']).area})
              #data.append({'field1':s1[field1], 'field2': s2[field2]})
    df = gpd.GeoDataFrame(data,columns=['geometry', field1, 'area'])
    if write == 'no':
        return df
    else:
        if sfname == 'default':
            outsfname = os.path.dirname(shp1) + '/' + os.path.splitext(os.path.basename(shp1))[0] + '.shp'
            df.to_file(outsfname)
            # Read in template proj
            prjref = shp1.replace(".shp",".prj")
            f = open(prjref)
            fr = f.read()
            f = None
            # Replace shp extension with prj extension for output file
            prjname = outsfname.replace(".shp",".prj")
            prj = open(prjname, "w")
            prj.write(fr)
            prj.close()
        else:
            outsfname = sfname
            df.to_file(outsfname)
            # Read in template proj
            prjref = shp1.replace(".shp",".prj")
            f = open(prjref)
            fr = f.read()
            f = None
            # Replace shp extension with prj extension for output file
            prjname = outsfname.replace(".shp",".prj")
            prj = open(prjname, "w")
            prj.write(fr)
            prj.close()
        return df
