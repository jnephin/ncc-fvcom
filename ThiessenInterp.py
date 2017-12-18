###############################################################################
#
# Author:       Jessica Nephin
# Affiliation:  IOS, Fisheries and Oceans Canada (DFO)
# Group:        Marine Spatial Ecology & Analysis, Ecosystems Science Division
# Address:      9860 West Saanich Road, Sidney, British Columbia, V8L 4B2, Canada
# Contact:      e-mail: jessica.nephin@dfo-mpo.gc.ca | tel: 250.363.6564
# Project:      Ocean Modelling Data Layers
# Code name:    ThiessenInterp.py
# Date started: Jan 23, 2017
# Date edited:  Jan 24, 2017
#
# Overview:
# Interpolate ocean current metrics (spring conditions and seasonal difference)
# from derived data originally from Pramod's FVCOM model. Use coastline to prevent
# extrapolation past land barriers.
#
# Requirements:
# The following input files in 'Data/Derived' directory
#   Currents.shp
#   SalTemp.shp
#
###############################################################################


# Import arcpy module
import arcpy			# import the ESRI arcpy module
import os               # import os module

# move up one directory
#os.chdir('F:/Abiotic_data/Ocean_Models/Pramod_2017/Scripts')
os.chdir('..')

# Check out the ArcGIS Spatial Analyst extension license
arcpy.CheckOutExtension("Spatial")
# Overwrite output
arcpy.env.overwriteOutput = True

# Path to ocean layers
ocean =  "Boundary/NCC_Area_BoP.shp"


# Interpolate #

shp = ["Data/Derived/SalTemp.shp","Data/Derived/Currents.shp"]
for s in shp:

	# Get Environmental variable name
	desc = arcpy.Describe(s)
	evname = desc.baseName 

	# Local variables:
	fc = "Data/Interp/"+evname+".gdb/"+evname+"_fc"
	Thiessen = "Data/Interp/"+evname+".gdb/"+evname+"_Thiessen"
	Thiessen_ocean = "Data/Interp/"+evname+".gdb/"+evname+"_Thiessen_ocean"
	Thiessen_single = "Data/Interp/"+evname+".gdb/"+evname+"_Thiessen_single"
	Thiessen_clipped = "Data/Interp/"+evname+".gdb/"+evname+"_Interp"

	# 1)
	# Create geodatabase file
	arcpy.CreateFileGDB_management ("Data/Interp", evname)

	# Add env. points to geodatabase
	arcpy.CopyFeatures_management (s, fc)

	# 2)
	# Create Thiessen Polygons
	arcpy.CreateThiessenPolygons_analysis (fc, Thiessen, "ALL")

	# 4)
	# Overlay with ocean layer
	arcpy.Intersect_analysis ([Thiessen, ocean], Thiessen_ocean)
	arcpy.MultipartToSinglepart_management(Thiessen_ocean, Thiessen_single)

	# 5)
	# Select Layer By Location
	# Remove polygons over 50 m from original points
	arcpy.MakeFeatureLayer_management (Thiessen_single, "tmp")
	arcpy.SelectLayerByLocation_management ("tmp", "WITHIN_A_DISTANCE", fc, "50 Meters", "NEW_SELECTION")
	arcpy.CopyFeatures_management ("tmp", Thiessen_clipped)
