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

# set workspace
os.chdir("E:/AbioticData/Ocean_Models/Pramod_NSB_FVCOM_2017")

# Set local variables
infeature = os.getcwd()+"/Data/Derived/Currents.shp"
zField = "spr_MnSp"
inBarrierFeature =  os.getcwd()+"/Boundary/CentralCoast.shp"
cellSize =  20.0
smoothing = 1

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")

# Execute Spline with Barriers
outSB = SplineWithBarriers(infeature, zField, inBarrierFeature, cellSize, smoothing)

# Use splitext to set the output name
outfilename = os.getcwd()+"/Data/Interp/CC_rmsBottomSpeed_20m.tif"

# Save the output
outSB.save(outfilename)
