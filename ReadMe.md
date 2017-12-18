Contact
--------
Author:       Jessica Nephin    
Affiliation:  IOS, Fisheries and Oceans Canada (DFO)    
Group:        Marine Spatial Ecology & Analysis, Ecosystems Science Division    
Address:      9860 West Saanich Road, Sidney, British Columbia, V8L 4B2, Canada    
Contact:      e-mail: jessica.nephin@dfo-mpo.gc.ca | tel: 250.363.6564    


Overview
--------
The bottom salinity/temperature and current data were sourced from DFO FVCOM circulation
model of the North Central Coast Area. For details on the model see the draft Techincal Report
in the 'Documentation' directory. The raw data originated from the stratified (baroclinic) model.
The model output was time-averaged over 29 days to get the conditions over a spring-neap tidal
cycle.


Methods
-------
Current velocities (u and v) can be used to represent different components of the currents
depending on the method of averaging. Root mean square of hourly velocities over a
29 day period was used to resolve tidal current speed. To represent the general circulation
component, speed was calculated from the mean of hourly velocities over a 29 day period.


Data
----
The sal/temp and current data was processed using 'ProcessModelData.R' to derive environmental
predictors for species distribution modelling. See Data/Derived directory for the processed data.
The derived sal/temp and current data was then interpolated using a Thiessen polygon or Spline
with barriers method (see ThiessenInterp.py and SplineBarriers.py). See 'Data/Interp' for the
interpolated layer. The interpolated layer located in the Data/Interp directory.

See Documentation/SalTemp_Summary.txt and Documentation/Currents_Summary.txt for a description
of each variable, including their range and units.


Caveats
-------
Interpolated data extends past the model domain in several inlets, rivers and estuaries.
The interpolated data in those regions should be used with caution.
