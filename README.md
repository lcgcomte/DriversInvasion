# DriversInvasion
R Scripts supporting the manuscript 'Parsing out the drivers of riverine fish invasion in the U.S.'

1.Script_Patterns_Invasion.R:
-Estimate life-history affinities
-Calculate ln-ratios density & life-history proportions
-Compute functional overlap
-Assess regional trends
-Null model random colonizations

2.Script_Patterns_HydroAlteration.R:
-Download flow time series from USGS website
-Prepare & clean flow time series
-Interpolate missing values (ARIMA)
-Backtransform & Change units of flow time series
-Smooth the flow time series using a 7-days running average
-Run wavelet analysis
-Computation flow metrics 

3.Script_Drivers_Invasion.R:
-Prepare datasets for models
-Run MCMCglmm models

--> associated with files:

#Occurrence database at the HUC8 scale (clean with names harmonized as described in Material & Methods):
"database_Fish.txt"

#Trait database (clean with imputed values as described in Material & Methods):
"database_traits_imputed.txt"

#File with characteristics of selected gages  [see SI Appendix, Table S2 for details on the variables]:
"Characteristics_selected_gages.txt"

#File with characteristics of sub-watersheds (=HUC8) [NHDPlus, available at http://nhd.usgs.gov/wbd.html]:
"WBD_HUC8.txt"

#Shapefile of local catchments (=HUC12) [NHDPlus, available at http://nhd.usgs.gov/wbd.html]:
"WDB_HUC12.shp"   

#File recreational freshwater fishing demand per HUC12 [available at https://www.epa.gov/enviroatlas/enviroatlas-data]:
"FreshwaterFishing_RecreationDemand.csv"





