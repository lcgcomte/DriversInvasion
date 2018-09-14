rm(list=ls(all=TRUE))
#///////////////////////////////////////////////////
#Download flow time series from USGS website
#///////////////////////////////////////////////////
require(dataRetrieval);

#File with characteristics of selected gages [see SI Appendix, Table S2 for details on variables]
gages = read.table("Characteristics_selected_gages.txt",h = T)
gages$Gage_ID = zeroPad(gages$Gage_ID,8)  #pad gage ID with zeros

for(i in 1:nrow(gages)){
siteNumber <-gages$Gage_ID[i]
parameterCd <- "00060" # discharge
rawDailyData <-readNWISdv(siteNumber, parameterCd, "1987-01-01","2016-12-31")
write.table(rawDailyData,paste("rawDailyData",siteNumber,".txt",sep=""),sep="\t")
print(i)
}

rm(list=ls(all=TRUE))
#//////////////////////////////////// //////////////////////
#Prepare & clean flow time series
#/////////////////////////////////////////////////////////
require(dataRetrieval);

NL = levels(factor(seq(as.Date("1987-01-01"), as.Date("2016-12-31"), by="days")))
gages = read.table("Characteristics_selected_gages.txt",h = T)
gages$Gage_ID = zeroPad(gages$Gage_ID,8)  #pad gage ID with zeros

Time_series = matrix(NA,length(NL),nrow(gages)) #create empty matrix
rownames(Time_series) = NL
colnames(Time_series) = gages$Gage_ID

for(i in 1:nrow(gages)){  
siteNumber <-gages$Gage_ID[i]
rawDailyData = try(read.table(paste("rawDailyData",siteNumber,".txt",sep=""),sep="\t"),silent=T)  

#identify column(s) with flow values
ind = grep("00060_00003$",colnames(rawDailyData))

#fix common problems
for(iii in 1:length(ind)){
rawDailyData[,ind[iii]][which(rawDailyData[,ind[iii]] < 0)] = NA  #put NA when the water flowed backward
rawDailyData[,ind[iii]][which(rawDailyData[,ind[iii]] == -999999)] = NA  #put NA when the water was icy
}
#merge flow columns  (if more than one)
if(length(ind) > 1){   
rawDailyData$X_00060_00003 = ifelse(is.na(rawDailyData[,ind[1]]) == F,rawDailyData[,ind[1]],rawDailyData[,ind[2]])
rawDailyData$X_00060_00003_cd = ifelse(is.na(rawDailyData[,ind[1]]) == F,as.character(rawDailyData[,ind[1]+1]),as.character(rawDailyData[,ind[2]+1]))
}
#remove duplicated days
dup = which(duplicated(rawDailyData$Date))
if(length(dup) >0){
rawDailyData = rawDailyData[-dup,]
}
#put daily values into empty matrix
fil = rawDailyData$X_00060_00003[match(rownames(Time_series),rawDailyData$Date)]
Time_series[,i] = fil

print(i)
}

write.table(Time_series,"All_time_series.txt",sep="\t")

rm(list=ls(all=TRUE))
#//////////////////////////////////////////////////////////////////////
#Interpolate missing values
#//////////////////////////////////////////////////////////////////////
require(forecast); 

Time_series = read.table("All_time_series.txt",h=T)
date_s = as.Date(rownames(Time_series))
Year = format(date_s,"%Y")
Jday = as.numeric(format(date_s, "%j"))

Time_series_inter = matrix(NA,nrow(Time_series),ncol(Time_series)) #create empty matrix

for(i in 1:ncol(Time_series)){ 

dat_imp = data.frame(y = log(Time_series[,i]+0.01,10),Year,Jday)      
y = ts(dat_imp$y)

#put raw values into empty matrix
Time_series_inter[,i] = dat_imp$y

#identify missing values
if(length(which(is.na(y))) > 0){

#fit ARIMA (Autoregressive Integrated Moving Average) model
fit = auto.arima(y,stepwise=FALSE, approximation=FALSE, D = 1) # D=1 for seasonal model

#apply Kalman filter
kr <- KalmanSmooth(y, fit$model)

#identify missing values to impute
id.na <- which(is.na(Time_series_inter[,i]))

#remove missing values at the begining and end of time series
aa = split(id.na, cumsum(c(1, diff(id.na) != 1)))
last = length(y)
first = 1
is.first = sapply(aa,function(x) length(which(x %in% first)))
is.last = sapply(aa,function(x) length(which(x %in% last)))
if(sum(is.first) > 0 | sum(is.last) > 0){
aa = aa[-c(which(is.last == 1),which(is.first == 1))]
id.na=unlist(aa)
}

#replace missing values by imputed values into matrix
for (ii in id.na){
Time_series_inter[ii,i] <- fit$model$Z %*% kr$smooth[ii,]
}
}
print(i)
}

colnames(Time_series_inter) = colnames(Time_series)
rownames(Time_series_inter) = rownames(Time_series)

write.table(Time_series_inter,"All_time_series_interpolated.txt",sep="\t")

rm(list=ls(all=TRUE))
#//////////////////////////////////////////////////////////////////////////
#Backtransform & Change units of flow time series
#//////////////////////////////////////////////////////////////////////////
Time_series_inter = read.table("All_time_series_interpolated.txt",h=T)
#back log-transform
Time_series_inter = apply(Time_series_inter,2,function(x) (10^x)-0.01)
#change potential values < 0 to 0 [due to imputation]
Time_series_inter = ifelse(Time_series_inter < 0,0,Time_series_inter) 
#transform from ft3/s to m3/s
Time_series_inter = apply(Time_series_inter,2,function(x) x * 0.0283168466)
write.table(Time_series_inter,"All_time_series_interpolated_detransformed_m3s.txt",sep="\t")

rm(list=ls(all=TRUE))
#///////////////////////////////////////////////////////////////////////////
#Smooth the flow time series using a 7-days running average 
#///////////////////////////////////////////////////////////////////////////
require(zoo);

Time_series_inter = read.table("All_time_series_interpolated_detransformed_m3s.txt",h = T)
Date = rownames(Time_series_inter)

rmean <- rollapplyr(Time_series_inter, 7, mean, na.rm = TRUE, fill = NA,align="center",by.column = TRUE)
rmean = data.frame(rmean)
rownames(rmean) = Date
write.table(rmean,"All_time_series_smoothed_m3s.txt",sep="\t")

rm(list=ls(all=TRUE))
#//////////////////////////////////////////////////////////////////////////////////
#Run wavelet analysis
#//////////////////////////////////////////////////////////////////////////////////
require(WaveletComp);require(dplyr)

Time_series_inter = read.table("All_time_series_smoothed_m3s.txt",h = T)

for(i in 17:ncol(Time_series_inter)){

#prepare time series
rawDailyData = data.frame(date = as.Date(rownames(Time_series_inter)), flow=Time_series_inter[,i])
rawDailyData = na.omit(rawDailyData) #remove the na at the beggining and end of time series

#run wavelet analysis
my.wt = analyze.wavelet(rawDailyData, "flow",method="AR", #select AR (=red noise)
                          loess.span=0,  #select no detrending 
                          dt=1, dj=1/20, 
                          lowerPeriod=2,
                          make.pval=T, n.sim=2,date.format = "%Y-%m-%d") 

#save the wavelet analysis outputs
save(my.wt, file = paste("Wavelet_Output",colnames(Time_series_inter)[i],".Rdata",sep="")) 
}

rm(list=ls(all=TRUE))
#///////////////////////////////////////////////////////////////////////////////////
#Computation flow metrics based on wavelet analysis
#//////////////////////////////////////////////////////////////////////////////////
require(WaveletComp);

Time_series_inter = read.table("All_time_series_smoothed_m3s.txt",h = T)

files = list.files(path = ".",pattern="Wavelet_Output")

SNR = AVP = NULL
for(i in 1:length(files)){
Name = gsub(".Rdata","",files[i])
Name = gsub("Wavelet_Output","",Name)
ind = grep(Name,colnames(Time_series_inter))
rawDailyData = data.frame(date = as.Date(rownames(Time_series_inter)), flow=Time_series_inter[,ind])
rawDailyData = na.omit(rawDailyData) #can remove the na because only at the beggining and end of time series

#load wavelet analysis outputs
load(paste("./",files[i],sep=""))

#identify scales < 365 days
freq = 1/my.wt$Period[c(which(my.wt$Period < 365.25),max(which(my.wt$Period < 365.25))+1)]

#identify scales over 365 days
freq.12 = 1/my.wt$Period[c(max(which(my.wt$Period < 365.25)),max(which(my.wt$Period < 365.25))+1)]

#////////////////////////////////
#Predictability: signal-to-noise ratio (SNR)
#////////////////////////////////
#select amplitudes of the global power spetrum for scales < 365 days
Amp = my.wt$Power.avg[c(which(my.wt$Period < 365.25),max(which(my.wt$Period < 365.25))+1)]

#calculate signal as the root mean square of significant frequencies
Amp.sign = ifelse(my.wt$Power.avg.pval[c(which(my.wt$Period < 365.25),max(which(my.wt$Period < 365.25))+1)] < 0.05,Amp,NA)
RMS.sign = sqrt(mean(c(Amp.sign),na.rm=T)) 

#calculate noise as the root mean square of non-significant frequencies
Amp.Nsign = ifelse(my.wt$Power.avg.pval[c(which(my.wt$Period < 365.25),max(which(my.wt$Period < 365.25))+1)] > 0.05,Amp,NA)
RMS.noise = sqrt(mean(c(Amp.Nsign),na.rm=T))

#calculate SNR
SNR = c(SNR,20*log(RMS.sign/RMS.noise,10))

#////////////////////////////////
#Seasonality: annual variation power in daily discharge (AVP) 
#////////////////////////////////
#select significant amplitudes over scales of 365 days for the entire time-series
Amp.12 = my.wt$Power[c(max(which(my.wt$Period < 365.25)),max(which(my.wt$Period < 365.25))+1),]

#integrate the scale-averaged wavelet power curve
y = apply(Amp.12,2,function(x)((12.77*1/20)/0.776)*weighted.mean(x,w = (1/freq.12),na.rm=T))
x = rawDailyData$date
AVP = c(AVP, sum(diff(as.numeric(x)) * (head(y,-1)+tail(y,-1)))/2)

print(i)
}

#save results
col.names = gsub("Wavelet_Output","",files)
col.names = gsub(".Rdata","",col.names)
names(AVP) = col.names
names(SNR) = col.names

write.table(data.frame(SNR),"SNR_wavelets.txt",sep="\t")
write.table(data.frame(AVP),"AVP_wavelets.txt",sep="\t")

rm(list=ls(all=TRUE))
#///////////////////////////////////////////////////////////////////////////////////
#Computation flow metrics based on median and CV
#//////////////////////////////////////////////////////////////////////////////////
require(dataRetrieval);

Time_series_inter = read.table("All_time_series_smoothed_m3s.txt",h = T)
Date = as.Date(rownames(Time_series_inter))

gages = read.table("Characteristics_selected_gages.txt",h = T)
gages$Gage_ID = zeroPad(gages$Gage_ID,8)  #pad gage ID with zeros

#////////////////////////////////
#Magnitude: annual median daily discharge standardized by the area of the watershed (Median) 
#////////////////////////////////
YR = sapply(strsplit(as.character(Date),"-"),'[',1)
Median.flow = apply(Time_series_inter,2,tapply,YR,median,na.rm=T)
Median.flow = t(Median.flow)
Median = data.frame(Median=apply(Median.flow,1,mean,na.rm=T))
#standardize Median by subwatershed area 
Median = Median/gages$WsAreaSqKm
#log-transform
Median = log(Median+0.00001,10)
write.table(Median,"Median_annual_flow_standardized.txt",sep="\t")

#////////////////////////////////
#Variability: annual coefficient of variation (CV) 
#////////////////////////////////
CV.flow = apply(Time_series_inter,2,function(x) tapply(x,YR,sd,na.rm=T)/tapply(x,YR,mean,na.rm=T))
CV.flow = t(CV.flow)
CV.flow[is.na(CV.flow)==T] = 0
CV = data.frame(CV=apply(CV.flow,1,mean,na.rm=T))
write.table(CV,"CV_annual_flow.txt",sep="\t")

rm(list=ls(all=TRUE))
#//////////////////////////////////////////////////
# RF models to hindcast synthetic hydrograph variables
#//////////////////////////////////////////////////

#////////////////////////////////////
#Prepare datasets
#////////////////////////////////////
require(dataRetrieval);

gages = read.table("Characteristics_selected_gages.txt",h = T)
gages$Gage_ID = zeroPad(gages$Gage_ID,8)  #pad gage ID with zeros

Median = read.table("Median_annual_flow_standardized.txt",h = T)
CV = read.table("CV_annual_flow.txt",h = T)
AVP = read.table("AVP_wavelets.txt",h = T)
SNR = read.table("SNR_wavelets.txt",h = T)

Y = data.frame(Median,CV,AVP,SNR)
rownames(Y) = gsub("X","",rownames(Y))
rownames(Y) = zeroPad(rownames(Y),8)  #pad gage ID with zeros
Y = data.frame(Gage_ID = rownames(Y),Y)

#predictors for calibration dataset
X.calib = data.frame(gages[which(gages$Category == "Reference"),-c(2:6)])  

#response variables for calibration dataset
Y.calib = Y[match(X.calib$Gage_ID,rownames(Y)),]

#predictors for prediction dataset
X.pred = data.frame(gages[which(gages$Category == "Altered"),-c(2:6)])

#response for prediction dataset (observed flow regime at altered gages)
Y.pred = data.frame(gages[which(gages$Category == "Altered"),-c(2:6)])
Y.pred = Y[match(X.pred$Gage_ID,rownames(Y)),]
Y.pred$Gage_ID = X.pred$Gage_ID #if some values are missing

write.csv(X.calib,"Predictors_calibration.csv", row.names=FALSE)
write.csv(Y.calib,"Responses_calibration.csv", row.names=FALSE)
write.csv(X.pred,"Predictors_prediction.csv", row.names=FALSE)
write.csv(Y.pred,"Responses_prediction.csv", row.names=FALSE)

rm(list=ls(all=TRUE))
#////////////////////////////////////
#Run models
#////////////////////////////////////
require(dplyr);require(randomForest);require(foreach);require(VSURF);require(doParallel)
registerDoParallel(cores=2)

response_data <- read.csv("Responses_calibration.csv")
predictor_data <- read.csv("Predictors_calibration.csv")
prediction_data <- read.csv("Predictors_prediction.csv")

for(i in 1:(ncol(response_data)-1)) {
  response_var <- names(response_data[i+1]) # name of response variable
  train <- response_data[c(1,(i+1))] # select Gage_ID and response variable data
  train <- left_join(train, predictor_data, by = "Gage_ID") # compile model training dataset (response variable and predictors)
  train <- na.omit(train) # remove observations with "NA" values
  sitelist <- as.character(unique(train$Gage_ID)) # list of sites

  #variable selection on full dataset
  RF = VSURF::VSURF(x = train[,3:ncol(train)], y=train[,2], ntree=2000,parallel = T)  #to select variables
  vari = RF$varselect.pred + 2 #+2 to take into account the fact that X start at the 3rd column

## Model training and performance assessment, based on 100 independent model iterations
for(indx in 1:100){

    #////////////////////////////////////////////////////////////////////
    #Cross-validation: select random sample of the sites per ecoregion respecting hydrologic unit prevalence
    #////////////////////////////////////////////////////////////////////
    tmp.idx = NULL
    for(kk in levels(train$ECOREG)){
    sitelist.eco = as.character(unique(train$Gage_ID[train$ECOREG == kk]))
    tmp.idx.eco <- as.character(sample(sitelist.eco,floor(0.8*length(sitelist.eco)),replace=FALSE)) # select random sample of sites (80%)
    tmp.idx = c(tmp.idx,tmp.idx.eco)
     }
    train.c <- train[train$Gage_ID %in% tmp.idx,] # calibration training set, excluding 20% of sites
    train.v <- train[!train$Gage_ID %in% tmp.idx,] # validation dataset, based on 20% of sites
    rm(tmp.idx)

    # train with calibration data and generate predictions at validation sites
    rf <- randomForest::randomForest(y = train.c[,2], x=train.c[,vari], importance = T, ntree=1000)  #with only a selection of variables
  
    ValPred <- predict(rf,train.v)

    # calculate model performance metrics
    validation <- data.frame(ValPred, ValObs = train.v[,2])
    colnames(validation) <- c("ValPred","ValObs")
    validation$resids <- validation$ValObs - validation$ValPred
    validation <- cbind(train.v[,1], validation) # add Gage_ID field

    results <- data.frame(Var = response_var,
                      RMSE = sqrt(mean((validation$ValObs - validation$ValPred)^2)), 
                      RSR = sqrt(mean((validation$ValObs - validation$ValPred)^2)) / sd(validation$ValObs), 
                      rsquared = cor(validation$ValPred,validation$ValObs)^2, 
                      pbias = ((sum(validation$ValObs - validation$ValPred))*100)/sum(validation$ValObs),
                      NSE = 1 - sum((validation$ValObs - validation$ValPred)^2) / sum((validation$ValObs - mean(validation$ValObs))^2),
                      MnOE = mean(validation$ValObs/validation$ValPred),
                      MedOE = median(validation$ValObs/validation$ValPred),
                      SDOE = sd(validation$ValObs/validation$ValPred))
    filename <- paste(response_var,"_", indx, "_model_perf.csv", sep="")
    write.csv(results, filename, row.names = F)
    
    preds<-predict(rf, prediction_data[,-1], type="response", predict.all = F) # save all model-tree predictions
    pred_results <- data.frame(Gage_ID = prediction_data[,1], preds)
    filename2 <- paste(response_var, "_",indx,"_predictions.csv", sep = "")
    write.csv(pred_results, filename2, row.names = F)

  }
}

rm(list=ls(all=TRUE))
#//////////////////////////////////////////
#Assess predictive performance
#//////////////////////////////////////////
response_data <- read.csv("Responses_calibration.csv")

for(i in 1:(ncol(response_data)-1)) {
response_var <- names(response_data)[i+1] # name of response variable
aa = list.files(".",paste("^",response_var,".*perf.csv$",sep=""))

perf_var = NULL
   for(j in 1:length(aa)){
   f = read.csv(aa[j])
   perf_var = rbind(perf_var,f[,-1])
   }
write.table(perf_var,paste0("Performance_models_",response_var,".txt"),sep="\t")
}

rm(list=ls(all=TRUE))
#//////////////////////////////////////////
#Aggregate predictions
#//////////////////////////////////////////
response_data <- read.csv("Responses_prediction.csv")

mean_pred = NULL
for(i in 1:(ncol(response_data)-1)) {
response_var <- names(response_data[i+1]) # name of response variable
aa = list.files(".",paste("^",response_var,".*predictions.csv$",sep=""))

pred_var = NULL
   for(j in 1:length(aa)){
   f = read.csv(aa[j])
   pred_var = cbind(pred_var,f[,-1])
   }
rownames(pred_var) = response_data$Gage_ID
write.table(data.frame(response_var=pred_var),paste0("Predictions_models_",response_var,".txt"),sep="\t")

#calculate mean prediction
mean_pred = cbind(mean_pred,apply(pred_var,1,mean))
}

colnames(mean_pred) = names(response_data)[-1] 
write.table(mean_pred,"Predictions_models_mean.txt",sep="\t")

rm(list=ls(all=TRUE))
#////////////////////////////////////////////
#Compute flow alterations [ln(obs/exp)]
#////////////////////////////////////////////
response_data <- read.csv("Responses_prediction.csv")
mean_pred = NULL
for(i in 1:(ncol(response_data)-1)) {
response_var <- names(response_data[i+1]) # name of response variable

hydro = read.table(paste0("Predictions_models_",response_var,".txt"),h = T)

aa = min(response_data[,response_var],unlist(hydro),na.rm=T)
if(aa <= 0){
eps = abs(min(c(response_data[,response_var],unlist(hydro)),na.rm=T))+0.01
alter = log((response_data[,response_var]+eps)/(hydro+eps))#add epsilon to avoid negative values
}else{
alter = log(response_data[,response_var]/hydro) 
}

rownames(alter) = response_data$Gage_ID
write.table(alter,paste0("Alteration_models_",response_var,".txt"),sep="\t")

#calculate mean alterations
mean_pred = cbind(mean_pred,apply(alter,1,mean))
}
colnames(mean_pred) = names(response_data)[-1] 
write.table(mean_pred,"Alteration_models_mean.txt",sep="\t")

rm(list=ls(all=TRUE))
#////////////////////////////////////////////////////////////////////
#Aggregate at the HUC8 scale
#////////////////////////////////////////////////////////////////////
require(dataRetrieval);

gages = read.table("Characteristics_selected_gages.txt",h = T)
gages$Gage_ID = zeroPad(gages$Gage_ID,8)  #pad gage ID with zeros

alt = read.table("Alteration_models_mean.txt",h = T)
rownames(alt) = zeroPad(rownames(alt),8)  #pad gage ID with zeros

#Weighted average at the HUC8 scale using local catchment areas as weights 
alt_huc8 = NULL 
for(i in 1:(ncol(alt))){
response_var <- names(alt[i]) # name of response variable

av_pre.alt = data.frame(Z = alt[,i],HUC8 = gages$HUC8[match(rownames(alt),gages$Gage_ID)],Surface = gages$CatAreaSqKm[match(rownames(alt),gages$Gage_ID)])
alt_huc8 = cbind(alt_huc8,sapply(split(av_pre.alt, av_pre.alt$HUC8), function(z) weighted.mean(z$Z, z$Surface)))
}
rownames(alt_huc8) = zeroPad(rownames(alt_huc8),8)  #pad gage ID with zeros
colnames(alt_huc8) = paste0("Alter_",names(alt))

#save file
write.table(alt_huc8,"Alterations_huc8.txt",sep="\t")

rm(list=ls(all=TRUE))
#/////////////////////////////////////////////////////////////////////
#Assess regional trends (HUC2 hydrologic unit)
#/////////////////////////////////////////////////////////////////////
Alter = read.table("Alterations_huc8.txt",h = T)

#create variable HUC2
HUC2 = sapply(strsplit(rownames(Alter),"*"),function(x) paste0(x[1],x[2],sep=""))

#Linear model without a separate intercept for each hydrologic unit
#Median
m.median = lm(Alter$Alter_Median ~ HUC2-1)
anova(m.median)
summary(m.median)

#CV
m.cv = lm(Alter$Alter_CV ~ HUC2-1)
anova(m.cv)
summary(m.cv)

#AVP
m.avp = lm(Alter$Alter_AVP ~ HUC2-1)
anova(m.avp)
summary(m.avp)

#SNR
m.snr = lm(Alter$Alter_SNR ~ HUC2-1)
anova(m.snr)                                        
summary(m.snr)
