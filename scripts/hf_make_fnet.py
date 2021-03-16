# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 20:55:08 2017
This script makes a fishnet.
@author: pj276
"""
import sys, os
sys.path.append('/home/pj276/projects/ifl_corridors/code/functions/')
#sys.path.append('/home/pj276/projects/ifl_corridors/code/parameter_files/')
#import p_121817 as p1
import hf_utilities as hfu

# Input directory
#indir = p1.infnetdir #'/gfc1_4_sa_ss'
#indir = '/scratch/pj276/ifl_corridors/fishnets'
# Output directory (edge matched tiles)
#odir = p1.outfnetdir #'/gfc1_4_sa_ss_emt'
#outdir = '/scratch/pj276/ifl_corridors/fishnets'
#outdir = '/scratch/pj276/undp_connectivity/fishnets'

# Fishnet file name
#fnet = p1.fnet
#fnet = '/projects/above_gedi/pjantz/ifl_corridors_data/fishnets/fnet_as_ss.shp'
#fnet = '/projects/above_gedi/pjantz/ifl_corridors_data/fishnets/fnet_sa_ss.shp'
#fnet = '/projects/above_gedi/pjantz/ifl_corridors_data/fishnets/fnet_sa_250k_ss.shp'
#fnet = '/projects/above_gedi/pjantz/ifl_corridors_data/fishnets/fnet_af_ss.shp'
#fnet = '/projects/above_gedi/pjantz/ifl_corridors_data/fishnets/fnet_sa_270k_ss.shp'
#fnet = '/projects/above_gedi/pjantz/ifl_corridors_data/fishnets/fnet_af_270k_ss.shp'
#fnet = '/scratch/pj276/fishnets/fnet_as_270k_ss.shp'
#fnet = '/scratch/pj276/fishnets/fnet_as_ss_v2.shp'

#fnet = '/scratch/pj276/undp_connectivity/fishnets/fnet_sa_ss_mspa.shp'
fnet = '/scratch/pj276/undp_connectivity/fishnets/fnet_af_ss_mspa.shp'
#fnet = '/scratch/pj276/undp_connectivity/fishnets/fnet_as_ss_mspa.shp'

# Reference raster for projection
#src_filename = p1.fnetref
#src_filename = '/projects/above_gedi/pjantz/ifl_corridors_data/gfc/gfc1_4_sa_ss/Hansen_GFC-2016-v1.4_lossyear_30N_090W_ss.tif'
# projection file name
#prjfile = '/projects/above_gedi/pjantz/ifl_corridors_data/prjfiles/ss_sa.prj'
#prjfile = '/projects/above_gedi/pjantz/ifl_corridors_data/prjfiles/ss_af.prj'
#prjfile = '/scratch/pj276/prjfiles/ss_as.prj'

#prjfile = '/scratch/pj276/ifl_corridors/prjfiles/ss_sa.prj'
prjfile = '/scratch/pj276/ifl_corridors/prjfiles/ss_af.prj'
#prjfile = '/scratch/pj276/ifl_corridors/prjfiles/ss_as.prj'

# If fishnet doesn't exist, create it
if not os.path.exists(fnet):
    # Upper left for South America and approximate xmax, ymin
    # This is for gfc15
#    xmin = -5755748.020845168270171
#    xmax = 3662811.9791548317
#    ymin = -3320236.6020596176
#    ymax = 3320293.397940382361412
    
    # This for gfc14
#    xmin = -3837548.020849999971688
#    xmax = 3162381.97915
#    ymin = -3679636.60206
#    ymax = 3320293.397940000053495

    # Upper left for Africa (from 10N_020W and 20N_030E)
    # This is for gfc15
#    xmin = -5009557.085697310976684
#    xmax = 3837562.914302689
#    ymin = -3320255.166765628
#    ymax = 2212524.833234372083098

    # This is for gfc14    
#    xmin = -5009557.085697310976684
#    xmax = 3722842.914302689
#    ymin = -3381425.166765628
#    ymax = 2212524.833234372083098

    # Asia
    # This is for gfc15
#    xmin = -5009551.383066884241998
#    xmax = 7126738.616933116
#    ymin = -3320263.745828367
#    ymax = 3320266.254171633161604
    
    # This is for gfc14
#    xmin = -4933771.383066884241998
#    xmax = 6946228.616933116
#    ymin = -3429883.745828367
#    ymax = 3320116.254171633161604

    # This is for mspa sa
#    xmin = -1860009.149951
#    xmax = 3362960.850049
#    ymin = -2276594.063266
#    ymax = 1392375.936734

    # This is for mspa sa expanded extent
#    xmin = -4948732.39861632
#    xmax = 4051177.6013836795
#    ymin = -3737348.16062553
#    ymax = 3262581.83937447

    # This is for mspa af
#    xmin = -2300332.542977
#    xmax = 2579047.457023
#    ymin = -896461.876964
#    ymax = 759958.123036

    # This is for mspa af expanded extent
    xmin = -4312435.51884177
    xmax = 4687474.48115823
    ymin = -3777923.0674729003
    ymax = 2222016.9325270997
        
    # This is for mspa as
#    xmin = -2747151.441167
#    xmax = 4277798.558833
#    ymin = -1278762.233736
#    ymax = 2864357.766264

        # This is for mspa as expanded extent
#    xmin = -4766943.90634788
#    xmax = 7232936.09365212
#    ymin = -3497347.87754044
#    ymax = 3502582.12245956
    
    # Cell height (cell dimensions may need to be a multiple of 30m)
    ch = 999990
    #ch = 500010 # for ~500km cells
    #ch = 270000 # for 270km cells (evenly divisible by 90)
    #ch = 249990 # for ~250km cells
    # Cell width
    cw = 999990
    #cw = 500010
    #cw = 270000 # for 270km cells
    #cw = 249990 # for ~250km cells
    # Get rows and columns by dividing the full extent of the fishnet in each dimension,
    # dividing by the grid cell dimension and adding 1 for slop. The units don't matter
    # just so long as they're consistent.
    # Or, can provide explicitly
    #cols = int(7000/500)+1
    #rows = int(7000/500)+1
    #cols = 10
    #rows = 7
    #cols = 9
    #rows = 6
    #cols = 13
    #rows = 7
    # sa mspa expanded
#    cols = 9
#    rows = 7
# af mspa expanded
    cols = 9
    rows = 6
# as mspa expanded
#    cols = 12
#    rows = 7
    #cols = 44
    #rows = 25
    #cols = 33
    #rows = 21
    #cols = 11
    #rows = 6
    #cols = 26 # for 270km cells
    #rows = 26 # for 270km cells
    #cols = 27
    #rows = 27
    #cols = 9
    #rows = 6
    # Get xmax and ymin
    #xmax = (math.ceil((axm-xmin)/cw)*cw)+xmin
    #ymin = ymax-(math.ceil((ymax-aym)/ch)*ch)
    # Create fishnet
    hfu.fishnet(fnet,xmin,ymax,rows,cols,str(ch),str(cw))
    #hfu.fishnet(fnet,xmin,xmax,ymin,ymax,str(ch), str(cw))
    # Get destination projection info
    # Get destination projection info
    f = open(prjfile)
    fr = f.read() # Read wkt to variable
    f = None # Close file
    # Create projection file name
    prjname = os.path.dirname(fnet) + os.sep + os.path.basename(fnet).replace("shp","prj")
    prj = open(prjname, "w")
    prj.write(fr)
    prj.close()
