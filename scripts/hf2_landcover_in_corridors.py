# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 15:15:49 2019

@author: pj276
"""

import arcpy, os
# Import spatial analyst
from arcpy.sa import *
# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")

# Corridor data year (should be 2012 but I messed up the naming)
# Don't change this without also changing the year for the land cover data
year = '2013'

# Africa
arcpy.env.extent = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_af_" + year + "_ifl_mos.tif"
arcpy.env.snapRaster = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_af_" + year + "_ifl_mos.tif"
arcpy.ProjectRaster_management("C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992_2015-v2.0.7.tif\Band_21",
                               "C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_af_ss.tif",
                               "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_af_" + year + "_ifl_mos.tif", "NEAREST", "90",
                               "#", "#", "#")
                               
# Asia
arcpy.env.extent = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_as_" + year + "_ifl_mos.tif"
arcpy.env.snapRaster = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_as_" + year + "_ifl_mos.tif"
arcpy.ProjectRaster_management("C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992_2015-v2.0.7.tif\Band_21",
                               "C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_as_ss.tif",
                               "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_as_" + year + "_ifl_mos.tif", "NEAREST", "90",
                               "#", "#", "#")

# South America
arcpy.env.extent = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_sa_" + year + "_ifl_mos.tif"
arcpy.env.snapRaster = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_sa_" + year + "_ifl_mos.tif"
arcpy.ProjectRaster_management("C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-300m-P1Y-1992_2015-v2.0.7.tif\Band_21",
                               "C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_sa_ss.tif",
                               "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_sa_" + year + "_ifl_mos.tif", "NEAREST", "90",
                               "#", "#", "#")

#---------------------------------
# Get land cover in corridor areas
# Africa
arcpy.env.extent = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_af_" + year + "_ifl_mos.tif"
arcpy.env.snapRaster = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_af_" + year + "_ifl_mos.tif"
outsa = Con(Raster("G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_af_" + year + "_ifl_mos.tif") > 0,
              "C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_af_ss.tif")
outsa.save("C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_af_ss_cmask.tif")

# Asia
arcpy.env.extent = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_as_" + year + "_ifl_mos.tif"
arcpy.env.snapRaster = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_as_" + year + "_ifl_mos.tif"
outsa = Con(Raster("G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_as_" + year + "_ifl_mos.tif") > 0,
              "C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_as_ss.tif")
outsa.save("C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_as_ss_cmask.tif")

# South America
arcpy.env.extent = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_sa_" + year + "_ifl_mos.tif"
arcpy.env.snapRaster = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_sa_" + year + "_ifl_mos.tif"
outsa = Con(Raster("G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_sa_" + year + "_ifl_mos.tif") > 0,
              "C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_sa_ss.tif")
outsa.save("C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_sa_ss_cmask.tif")

#---------------------------------
# Mask out ifls
# Africa
arcpy.env.extent = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_af_" + year + "_ifl_mos.tif"
arcpy.env.snapRaster = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_af_" + year + "_ifl_mos.tif"
ifl = "G:/My Drive/Projects/ifl_corridors/geodata/ifls/IFL_2000/ifl_2000_tropics_ssu_af.shp"
if not os.path.exists(ifl):
    arcpy.FeatureToRaster_conversion(ifl, "FID", "G:/My Drive/Projects/ifl_corridors/geodata/ifls/IFL_2000/ifl_2000_tropics_ssu_af_90m.tif", 90)
outsa = Con(IsNull(Raster("G:/My Drive/Projects/ifl_corridors/geodata/ifls/IFL_2000/ifl_2000_tropics_ssu_af_90m.tif")),
           "C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_af_ss_cmask.tif")
outsa.save("C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_af_ss_iflcmask.tif")

# Asia
arcpy.env.extent = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_as_" + year + "_ifl_mos.tif"
arcpy.env.snapRaster = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_as_" + year + "_ifl_mos.tif"
ifl = "G:/My Drive/Projects/ifl_corridors/geodata/ifls/IFL_2000/ifl_2000_tropics_ssu_as.shp"
if not os.path.exists(ifl):
    arcpy.FeatureToRaster_conversion(ifl, "FID", "G:/My Drive/Projects/ifl_corridors/geodata/ifls/IFL_2000/ifl_2000_tropics_ssu_as_90m.tif", 90)
outsa = Con(IsNull(Raster("G:/My Drive/Projects/ifl_corridors/geodata/ifls/IFL_2000/ifl_2000_tropics_ssu_as_90m.tif")),
           "C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_as_ss_cmask.tif")
outsa.save("C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_as_ss_iflcmask.tif")

# South America
arcpy.env.extent = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_sa_" + year + "_ifl_mos.tif"
arcpy.env.snapRaster = "G:/My Drive/Projects/ifl_corridors/rp5k130k25k/csum_rp5k130k90m_neg01_sa_" + year + "_ifl_mos.tif"
ifl = "G:/My Drive/Projects/ifl_corridors/geodata/ifls/IFL_2000/ifl_2000_tropics_ssu_sa.shp"
if not os.path.exists(ifl):
    arcpy.FeatureToRaster_conversion(ifl, "FID", "G:/My Drive/Projects/ifl_corridors/geodata/ifls/IFL_2000/ifl_2000_tropics_ssu_sa_90m.tif", 90)
outsa = Con(IsNull(Raster("G:/My Drive/Projects/ifl_corridors/geodata/ifls/IFL_2000/ifl_2000_tropics_ssu_sa_90m.tif")),
           "C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_sa_ss_cmask.tif")
outsa.save("C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_sa_ss_iflcmask.tif")

#---------------------------------
# Save raster attribute tables as csvs
# Africa
arcpy.TableToTable_conversion("C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_af_ss_iflcmask.tif",
                              "G:/My Drive/Projects/ifl_corridors/statsfiles/lc_stats", "esa_2012_lc_in_af_neg01_corridors_no_ifl.csv")
# Asia
arcpy.TableToTable_conversion("C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_as_ss_iflcmask.tif",
                              "G:/My Drive/Projects/ifl_corridors/statsfiles/lc_stats", "esa_2012_lc_in_as_neg01_corridors_no_ifl.csv")
# South America
arcpy.TableToTable_conversion("C:/data/landcover/esa_cci_global_300/ESACCI-LC-L4-LCCS-Map-90m-P1Y-2012-v2.0.7_sa_ss_iflcmask.tif",
                              "G:/My Drive/Projects/ifl_corridors/statsfiles/lc_stats", "esa_2012_lc_in_sa_neg01_corridors_no_ifl.csv")



