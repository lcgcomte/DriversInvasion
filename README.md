# DriversInvasion
R Scripts supporting the manuscript 'Parsing out the drivers of riverine fish invasion in the U.S.'

1.Script_Patterns_Invasion.R
-Estimate life-history affinities
-Calculate ln-ratios density & life-history proportions
-Compute functional overlap
-Assess regional trends
-Null model random colonizations

2.Script_Patterns_HydroAlteration.R
-Download flow time series from USGS website
-Prepare & clean flow time series
-Interpolate missing values (ARIMA)
-Backtransform & Change units of flow time series
-Smooth the flow time series using a 7-days running average
-Run wavelet analysis
-Computation flow metrics 

3.Script_Drivers_Invasion.R
-Prepare datasets for models
-Run MCMCglmm models
