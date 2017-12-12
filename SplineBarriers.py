# Name: SplineBarriers.py
# Description: Interpolate a series of point features onto a
#    rectangular raster using a barrier, using a
#    minimum curvature spline technique.
# Requirements: Spatial Analyst Extension and Java Runtime
# Author: Jessica Nephin

# Import system modules
import os
import arcpy
from arcpy import env
from arcpy.sa import *

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")

# set workspace
os.chdir("..")
#os.chdir("E:/AbioticData/Ocean_Models/Pramod_NSB_FVCOM_2017")

# Set local variables
infeature = os.getcwd()+"/Data/Derived/Currents.shp"

# Create a list of fields using the ListFields function
fields = arcpy.ListFields(infeature,"","Double")

# loop through numeric fields
for f in fields:

    # inputs
    zField = f.name
    inBarrierFeature =  os.getcwd()+"/Boundary/NCC_Area_20mbuffer_15ksponge.shp"
    cellSize =  20.0
    smoothing = 1

    # Execute Spline with Barriers
    outSB = SplineWithBarriers(infeature, zField, inBarrierFeature, cellSize, smoothing)

    # Use splitext to set the output name
    outfilename = os.getcwd()+"/Data/Interp/Rasters/FVCOM_"+zField+"_20m.tif"

    # Save the output
    outSB.save(outfilename)
