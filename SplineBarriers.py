# Name: SplineBarriers.py
# Description: Interpolate a series of point features onto a
#    rectangular raster using a barrier, using a
#    minimum curvature spline technique.
# Requirements: Spatial Analyst Extension and Java Runtime
# Author: Jessica Nephin

# Run from parent directory

# Import system modules
import os
import arcpy
from arcpy import env
from arcpy.sa import *

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")

# set workspace
#os.chdir("E:/AbioticData/Ocean_Models/Pramod_NSB_FVCOM_2017")

# Set local variables
infeature = os.getcwd()+"/Data/Derived/Currents.shp"

# Create a list of fields using the ListFields function
fields = arcpy.ListFields(infeature,"","Double")

# loop through numeric fields
for f in fields:

    # inputs
    zField = f.name
    inBarrierFeature =  os.getcwd()+"/Boundary/NCC_Coastline_Simplified_15kmsponge_50mbuffer.shp"
    inMaskData = os.getcwd()+"/Boundary/NCC_Area_15ksponge.shp"
    outfilename = os.getcwd()+"/Data/Interp/Rasters/FVCOM_"+zField+"_20m.tif"
    cellSize =  100.0
    smoothing = 1

    # Execute Spline with Barriers
    outSB = SplineWithBarriers(infeature, zField, inBarrierFeature, cellSize, smoothing)

    # Execute ExtractByMask
    outExtractByMask = ExtractByMask(outSB, inMaskData)

    # Save the output
    outExtractByMask.save(outfilename)
