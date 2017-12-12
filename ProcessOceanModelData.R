###############################################################################
#
# Author:       Jessica Nephin
# Affiliation:  IOS, Fisheries and Oceans Canada (DFO)
# Group:        Marine Spatial Ecology & Analysis, Ecosystems Science Division
# Address:      9860 West Saanich Road, Sidney, British Columbia, V8L 4B2, Canada
# Contact:      e-mail: jessica.nephin@dfo-mpo.gc.ca | tel: 250.363.6564
# Project:      Ocean Modelling Data Layers
# Code name:    ProcessOceanModelData.R
# Date started: Jan 23, 2017
# Date edited:  Jan 24, 2017
#
# Overview:
# Calculate ocean current metrics (spring conditions and seasonal difference)
# from raw data output from Pramod's FVCOM model. Project and clip points to
# model domain then export shapefiles. Export summary of data.
#
# Requirements:
# The following input files in 'Data/Raw' directory
#   element-locations.dat
#   fall-stratified-bottom_current.dat
#   fall-stratified-tauc.dat
#   fall-temp_sal.dat
#   spring-stratified-bottom_current.dat
#   spring-stratified-tauc.dat
#   spring-temp_sal.dat
#   summer-stratified-bottom_current.dat
#   summer-stratified-tauc.dat
#   summer-temp_sal.dat
#
###############################################################################


# move back to parent directory
setwd('..')

# packages
require(readr)
require(dplyr)
require(sp)
require(rgdal)
require(rgeos)
require(reshape2)


# load elements for current and tauc lat and long data
elm <- "Data/Raw/element-locations.dat"
elm <- read_table(file=elm, col_names = FALSE)


#################
#   Currents    #
#################

# load current and shear stress data
current_files <- list.files(path = "Data/Raw", pattern = "tauc|current", full.names = TRUE)
raw <- data.frame(ID = 1:209890)
for( i in current_files){
  tmp <- read_table(file=i)
  season <- sub(".*/","",i)
  season <- sub("-.*","",season)
  names(tmp) <- paste(season, names(tmp), sep = ".")
  raw <- data.frame(raw,tmp)
}

# add lat and lon
raw <- data.frame(raw, Longitude = elm$X1, Latitude = elm$X2)

# calculate and rename
cur <- raw %>%
  mutate(
    # circulation speed for each season
    spr.cur=sqrt(spring.umean^2 + spring.vmean^2),
    sum.cur=sqrt(summer.umean^2 + summer.vmean^2),
    fal.cur=sqrt(fall.umean^2 + fall.vmean^2)) %>%
  #select
  select( Longitude, Latitude,
          spr.cur, sum.cur, fal.cur, 
          spring.RMSspeed, summer.RMSspeed, fall.RMSspeed, 
          spring.stress, summer.stress, fall.stress) %>%
  rename( spr.rms = spring.RMSspeed, sum.rms = summer.RMSspeed, fal.rms = fall.RMSspeed,
          spr.stres = spring.stress, sum.stres = summer.stress, fal.stres = fall.stress) %>%
  data.frame( )

# plot and check for correlation
#plot(cur[,3:11])
#cor(cur[,3:11])



################
#   SalTemp    #
################

# load sal and temp data
saltemp_files <- list.files(path = "Data/Raw", pattern = "temp_sal", full.names = TRUE)
raw <- data.frame(ID = 1:118089)
for( i in saltemp_files){
  tmp <- read.table(file=i, header=F, sep=" ")
  season <- sub(".*/","",i)
  season <- sub("-.*","",season)
  names(tmp) <- c("Longitude", "Latitude", "Temp", "Sal")
  names(tmp) <- paste(season, names(tmp), sep = ".")
  raw <- data.frame(raw,tmp)
}

# calculate and rename
saltemp <- raw %>%
  rename( Longitude=fall.Longitude, Latitude=fall.Latitude ) %>%
  select( Longitude, Latitude,
          spring.Sal, summer.Sal, fall.Sal,
          spring.Temp, summer.Temp, fall.Temp ) %>%
  rename( spr.sal = spring.Sal, sum.sal = summer.Sal, fal.sal = fall.Sal,
          spr.temp = spring.Temp, sum.temp = summer.Temp, fal.temp = fall.Temp ) %>%
  data.frame( )


# plot and check for correlation
#plot(saltemp[,3:6])
#cor(saltemp[,3:6])




#################################
#    Project and Clip points    #
#################################

# create spatial point dataframes
coordinates(cur) <- ~Longitude+Latitude
coordinates(saltemp) <- ~Longitude+Latitude

#full bc albers proj
proj <- "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

# project to albers
proj4string(cur) <-CRS("+proj=longlat +datum=WGS84")
proj4string(saltemp) <-CRS("+proj=longlat +datum=WGS84")
cur <- spTransform(cur, proj)
saltemp <- spTransform(saltemp, proj)

# load model domain with 15km buffer for sponge layer
sponge <- readOGR(dsn = "Boundary", layer = "Edges_15km")
sponge <- spTransform(sponge, proj)

# Get the edges (sponge layer)
sponge <- readOGR( dsn="Boundary", layer="Edges_15km", verbose=FALSE )
# Set the projection
sponge <- spTransform(sponge, proj)

# Determine which points are within the sponge layer
inSponge.cur <- over(cur, sponge)
inSponge.cur <- is.na(inSponge.cur$BUFF_DIST)
inSponge.saltemp <- over(saltemp, sponge)
inSponge.saltemp <- is.na(inSponge.saltemp$BUFF_DIST)

# Remove points that are not within the sponge layer
cur.erase <- cur[inSponge.cur, ]
saltemp.erase <- saltemp[inSponge.saltemp, ]

#check
#plot(cur.erase)
#plot(sponge, add = T)
#plot(saltemp.erase)
#plot(sponge, add = T)

# write layers
writeOGR(cur.erase, dsn = "Data/Derived", layer = "Currents", 
         driver = "ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(saltemp.erase, dsn = "Data/Derived", layer = "SalTemp", 
         driver = "ESRI Shapefile", overwrite_layer = TRUE)




###################################
#      Sal/Temp Data Summary      #
###################################

# project to lat and lon
saltemp <- spTransform(saltemp.erase, CRS("+proj=longlat +datum=WGS84"))

# move back into a data.frame
saltemp <- as.data.frame(saltemp)

# summary
saltemp.melt <- melt(saltemp)
saltemp.table <- saltemp.melt %>%
  group_by(variable) %>%
  summarise(Min = round(min(value),2),
            Mean = round(mean(value),2),
            Max = round(max(value),2)) %>%
  as.data.frame()

# add units
saltemp.table$Units <- c("PSU",
                         "PSU",
                         "PSU",
                         "Celsius",
                         "Celsius",
                         "Celsius",
                         "Decimal degree",
                         "Decimal degree")

# add description
saltemp.table$Description <- c("Mean spring salinity",
                               "Mean summer salinity",
                               "Mean fall salinity",
                               "Mean spring temperature",
                               "Mean summer temperature",
                               "Mean fall temperature",
                               "Longitude (WGS 84)",
                               "Latitude (WGS 84)")


# set string length and alinement
saltemp.table$variable <- sprintf("%-10s", as.character(saltemp.table$variable))
saltemp.table$Min <- sprintf("%-8s", as.character(saltemp.table$Min))
saltemp.table$Mean <- sprintf("%-8s", as.character(saltemp.table$Mean))
saltemp.table$Max <- sprintf("%-8s", as.character(saltemp.table$Max))
saltemp.table$Units <- sprintf("%-15s", as.character(saltemp.table$Units))
saltemp.table$Description <- sprintf("%-70s", as.character(saltemp.table$Description))
names(saltemp.table) <- sprintf("%-8s", names(saltemp.table))
saltemp.table

# export summary
header <- "This is a summary of the processed bottom salinity and temperature data which can be found in 'Data/Derived' directory
\n"
notes <- "
Note: Created using  ProcessOceanModelData.R"
cat(header, file="Documentation/SalTemp_Summary.txt")
write.table(saltemp.table, file="Documentation/SalTemp_Summary.txt", append=TRUE,
            quote = FALSE, sep = "\t", row.name = FALSE)
cat(notes, file="Documentation/SalTemp_Summary.txt", append=TRUE)





###################################
#      Currents Data Summary      #
###################################

# project to lat and lon
cur <- spTransform(cur.erase, CRS("+proj=longlat +datum=WGS84"))

# move back into a data.frame
cur <- as.data.frame(cur)

# summary
cur.melt <- melt(cur)
cur.table <- cur.melt %>%
  group_by(variable) %>%
  summarise(Min = round(min(value),2),
            Mean = round(mean(value),2),
            Max = round(max(value),2)) %>%
  as.data.frame()

# add units
cur.table$Units <- c("m/s",
                     "m/s",
                     "m/s",
                     "m/s",
                     "m/s",
                     "m/s",
                     "log10(m2/s2)",
                     "log10(m2/s2)",
                     "log10(m2/s2)",
                     "Decimal degree",
                     "Decimal degree")

# add description
cur.table$Description <- c("Mean circulation spring bottom current speed",
                           "Mean circulation summer bottom current speed",
                           "Mean circulation fall bottom current speed",
                           "Root mean square spring bottom current speed",
                           "Root mean square summer bottom current speed",
                           "Root mean square fall bottom current speed",
                           "Mean spring shear stress",
                           "Mean summer shear stress",
                           "Mean fall shear stress",
                               "Longitude (WGS 84)",
                               "Latitude (WGS 84)")

# set string length and alinement
cur.table$variable <- sprintf("%-10s", as.character(cur.table$variable))
cur.table$Min <- sprintf("%-8s", as.character(cur.table$Min))
cur.table$Mean <- sprintf("%-8s", as.character(cur.table$Mean))
cur.table$Max <- sprintf("%-8s", as.character(cur.table$Max))
cur.table$Units <- sprintf("%-15s", as.character(cur.table$Units))
cur.table$Description <- sprintf("%-70s", as.character(cur.table$Description))
names(cur.table) <- sprintf("%-8s", names(cur.table))
cur.table

# export summary
header <- "This is a summary of the processed bottom current data which can be found in 'Data/Derived' directory
\n"
notes <- "
Note: Created using script ProcessOceanModelData.R"
cat(header, file="Documentation/Currents_Summary.txt")
write.table(cur.table, file="Documentation/Currents_Summary.txt", append=TRUE,
            quote = FALSE, sep = "\t", row.name = FALSE)
cat(notes, file="Documentation/Currents_Summary.txt", append=TRUE)




