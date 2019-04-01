rm(list=ls(all=TRUE))
#/////////////////////////////////////////////////////////////////////
#Prepare datasets for models
#/////////////////////////////////////////////////////////////////////
require(dataRetrieval);

gages = read.table("Characteristics_selected_gages.txt",h = T)
gages$Gage_ID = zeroPad(gages$Gage_ID,8)  #pad gage ID with zeros

#////////////////
#1-Niche oppurtunities arising from human-induced flow alteration (already HUC8 scale)
#////////////////
alt_huc8 = read.table("Alterations_huc8.txt",h=T)
rownames(alt_huc8) = zeroPad(rownames(alt_huc8),8)  #pad gage ID with zeros

#////////////////
#3-Propagule pressure
#////////////////
require(dams); require(rgdal)

#Import NHDPlus shapefile of local catchments (=HUC12) [available at http://nhd.usgs.gov/wbd.html]
WDB_HUC12 = readOGR(".","WDB_HUC12")   
crs_NAD83 = proj4string(WDB_HUC12) #save coordinate system

#Recreational freshwater fishing demand
#Import file per HUC12 [available at https://www.epa.gov/enviroatlas/enviroatlas-data]
Fishing = read.table("FreshwaterFishing_RecreationDemand.csv",sep=",",h=T)
Fishing$HUC_12 = zeroPad(Fishing$HUC_12,12)  #pad with zeros

#Weighted average at the HUC8 scale using local catchment areas as weights
Fishing = data.frame(Fishing,HUC8 = WDB_HUC12@data$HUC_8[match(Fishing$HUC_12,WDB_HUC12$HUC_12)],Surface = WDB_HUC12@data$AreaHUC12[match(Fishing$HUC_12,WDB_HUC12$HUC_12)])
Fishing_huc8 = sapply(split(Fishing, Fishing$HUC8), function(z) weighted.mean(z$FF_Demand, z$Surface))

#Densities of dams
#Select large dams [dams 50 feet or more in height, dams with a normal storage capacity of 5,000 acre-feet or more, and dams with a maximum storage capacity of 25,000 acre-feet or more]
data(nid_cleaned); dams = nid_cleaned
Height = apply(dams[,21:24],1,max,na.rm=T)
LD = dams[which(Height >= 50 & dams$NID_Storage >= 5000 | dams$NID_Storage >= 25000),]

#Identify the local catchment dams belong tot
LD = LD[-which(is.na(LD$Latitude)),]  #remove dams with unknown coordinates
coor_LD = SpatialPoints(LD[,6:7])
proj4string(coor_LD) <- crs_NAD83 #assign spatial coordinate system
ov_LD = over(coor_LD,WDB_HUC12)  #perform spatial join based on coordinates
NB_LD = table(ov_LD$HUC_12)

#Compute densities of large dams per local ctachment 
NB_LD_dens = NB_LD/WDB_HUC12$AreaHUC12[match(names(NB_LD),WDB_HUC12$HUC_12)]

#Weighted average at the HUC8 scale using local catchment areas as weights
dens_LD = data.frame(Z = NB_LD_dens,HUC12 = names(NB_LD_dens),HUC8 = sapply(strsplit(names(NB_LD_dens),"*"),function(x) paste0(x[1:8],collapse="")),Surface = WDB_HUC12$AreaHUC12[match(names(NB_LD_dens),WDB_HUC12$HUC_12)])  
dens_LD_huc8 = sapply(split(dens_LD, dens_LD$HUC8), function(z) weighted.mean(z$Z.Freq, z$Surface))

#//////////////////////////////////////
#Creating predictor file for models
#/////////////////////////////////////
dat.mod = data.frame(alt_huc8,Fishing = Fishing_huc8[match(rownames(alt_huc8),names(Fishing_huc8))],Dams = dens_LD_huc8[match(rownames(alt_huc8),names(dens_LD_huc8))])

#Save file
write.table(dat.mod,"Predictors_NAS_HUC8.txt",sep="\t")

rm(list=ls(all=TRUE));

#/////////////////////////////////////////////////////////////////////
#Run MCMCglmm models
#/////////////////////////////////////////////////////////////////////
require(MCMCglmm); require(coda)

#predictor variables
dat.mod = read.table("Predictors_NAS_HUC8.txt",h=T)
dat.mod$NB_LD_dens_huc8 = log(dat.mod$Dams + 0.0001)
dat.mod$Fishing = log(dat.mod$Fishing)

#responses variables
lnratio_prop=read.table("LnRatio_proportion_life-history_NAS_NA.txt",h=T)
lnratio_dens=read.table("LnRatio_densities_NAS_NA.txt",h=T)

#///////////////////////////////////////////////////////////////
#function to transform to z-scores
zscores = function(x) {Z = (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
return(Z)}
#////////////////////////////////////////////////////////////////

##/////////////////////////////
#Model ln-ratio densities (run 3 chains)
#/////////////////////////////
Y = lnratio_dens[,1] 
X = dat.mod[match(rownames(lnratio_dens),rownames(dat.mod)),]
DAT = data.frame(Y,X)
DAT = apply(DAT,2,zscores)  #transform to z-scores (both predictors and response)
DAT = data.frame(DAT)
DAT.dens = na.omit(DAT)

equa = "Y ~ Alter_SNR * Fishing + Alter_AVP * Fishing + Alter_CV * Fishing + Alter_Median * Fishing + Alter_SNR * Dams + Alter_AVP * Dams + Alter_CV * Dams + Alter_Median * Dams"

mm.ratio_tot1 = MCMCglmm(as.formula(equa),    
                    data=DAT.dens,nitt= 350000,thin=200,burnin= 60000, 
                    rcov = ~ units,     
                    verbose=FALSE,family=c("gaussian"))
                     
mm.ratio_tot2 = MCMCglmm(as.formula(equa),    
                    data=DAT.dens,nitt= 350000,thin=200,burnin= 60000, 
                    rcov = ~ units,     
                    verbose=FALSE,family=c("gaussian"))

mm.ratio_tot3 = MCMCglmm(as.formula(equa),    
                    data=DAT.dens,nitt= 350000,thin=200,burnin= 60000, 
                    rcov = ~ units,   
                    verbose=FALSE,family=c("gaussian"))

#//////////////////////
#Diagnostics
#///////////////////////
#Effective size
which(c(effectiveSize(mm.ratio_tot1$VCV),effectiveSize(mm.ratio_tot1$Sol))< 1000)#ok
which(c(effectiveSize(mm.ratio_tot2$VCV),effectiveSize(mm.ratio_tot2$Sol))< 1000) #ok
which(c(effectiveSize(mm.ratio_tot2$VCV),effectiveSize(mm.ratio_tot3$Sol))< 1000) #ok

#Autocorrelation
which(c(autocorr(mm.ratio_tot1$VCV)[2,,],autocorr(mm.ratio_tot1$Sol)[2,,]) > 0.1)#ok
which(c(autocorr(mm.ratio_tot2$VCV)[2,,],autocorr(mm.ratio_tot2$Sol)[2,,]) > 0.1) #ok
which(c(autocorr(mm.ratio_tot3$VCV)[2,,],autocorr(mm.ratio_tot3$Sol)[2,,]) > 0.1) #ok

#Heidel test
which(c(heidel.diag(mm.ratio_tot1$VCV)[,3],heidel.diag(mm.ratio_tot1$Sol)[,3]) < 0.05) #ok
which(c(heidel.diag(mm.ratio_tot2$VCV)[,3],heidel.diag(mm.ratio_tot2$Sol)[,3]) < 0.05) #ok
which(c(heidel.diag(mm.ratio_tot3$VCV)[,3],heidel.diag(mm.ratio_tot3$Sol)[,3]) < 0.05) #ok

#Mixing chains
m_tot_Sol <- mcmc.list(mm.ratio_tot1$Sol,mm.ratio_tot2$Sol,mm.ratio_tot3$Sol)
m_tot_VCV <- mcmc.list(mm.ratio_tot1$VCV,mm.ratio_tot2$VCV,mm.ratio_tot3$VCV)
aa.dens <- as.mcmc(do.call(rbind,m_tot_Sol))
bb.dens <- as.mcmc(do.call(rbind,m_tot_VCV))

#Gelman test
which(c(gelman.diag(m_tot_Sol,multivariate=F)$psrf[,2],gelman.diag(m_tot_VCV,multivariate=F)$psrf[,2]) >= 1.1) #ok

#Save output
write.table(aa.dens,"Coefficients_densities_MCMCglmm.txt",sep="\t")

##/////////////////////////////
#Model ln-ratio life-history proportions (run 3 chains)
#/////////////////////////////
Y = lnratio_prop 
X = dat.mod[match(rownames(lnratio_prop),rownames(dat.mod)),]
DAT = data.frame(Y,X)
DAT = apply(DAT,2,zscores)  #transform to z-scores (both predictors and response)
DAT = data.frame(DAT)
DAT.prop = na.omit(DAT)

equa = "cbind(equilibrium,opportunistic,periodic) ~ (trait - 1)+
trait:(Alter_SNR * Fishing + Alter_AVP * Fishing + Alter_CV * Fishing + Alter_Median * Fishing + Alter_SNR * Dams + Alter_AVP * Dams + Alter_CV * Dams + Alter_Median * Dams)"

mm.ratio.prop1 = MCMCglmm(as.formula(equa),  
                 data=DAT.prop,nitt= 350000,thin=200,burnin= 60000, 
                 rcov = ~ us(trait):units,
                 verbose=FALSE,family=c("gaussian","gaussian","gaussian"))

mm.ratio.prop2 = MCMCglmm(as.formula(equa),  
                 data=DAT.prop,nitt= 350000,thin=200,burnin= 60000, 
                 rcov = ~ us(trait):units,
                 verbose=FALSE,family=c("gaussian","gaussian","gaussian"))

mm.ratio.prop3 = MCMCglmm(as.formula(equa),  
                 data=DAT.prop,nitt= 350000,thin=200,burnin= 60000, 
                 rcov = ~ us(trait):units,
                 verbose=FALSE,family=c("gaussian","gaussian","gaussian"))

#//////////////////////
#Diagnostics
#///////////////////////
#Effective size
which(c(effectiveSize(mm.ratio.prop1$VCV),effectiveSize(mm.ratio.prop1$Sol))< 1000) #ok 
which(c(effectiveSize(mm.ratio.prop2$VCV),effectiveSize(mm.ratio.prop2$Sol))< 1000)##ok
which(c(effectiveSize(mm.ratio.prop3$VCV),effectiveSize(mm.ratio.prop3$Sol))< 1000)#ok

#Autocorrelation
which(c(autocorr(mm.ratio.prop1$VCV)[2,,],autocorr(mm.ratio.prop1$Sol)[2,,]) > 0.1) #ok
which(c(autocorr(mm.ratio.prop2$VCV)[2,,],autocorr(mm.ratio.prop2$Sol)[2,,]) > 0.1) #ok
which(c(autocorr(mm.ratio.prop3$VCV)[2,,],autocorr(mm.ratio.prop3$Sol)[2,,]) > 0.1) #ok

#Heidel test
which(c(heidel.diag(mm.ratio.prop1$VCV)[,3],heidel.diag(mm.ratio.prop1$Sol)[,3]) < 0.05) #ok
which(c(heidel.diag(mm.ratio.prop2$VCV)[,3],heidel.diag(mm.ratio.prop2$Sol)[,3]) < 0.05) #ok
which(c(heidel.diag(mm.ratio.prop3$VCV)[,3],heidel.diag(mm.ratio.prop3$Sol)[,3]) < 0.05) #ok

#Mixing chains
m_tot_Sol <- mcmc.list(mm.ratio.prop1$Sol,mm.ratio.prop2$Sol,mm.ratio.prop3$Sol)
m_tot_VCV <- mcmc.list(mm.ratio.prop1$VCV,mm.ratio.prop2$VCV,mm.ratio.prop3$VCV)
aa.prop <- as.mcmc(do.call(rbind,m_tot_Sol))
bb.prop <- as.mcmc(do.call(rbind,m_tot_VCV))

#Gelman test
which(c(gelman.diag(m_tot_Sol,multivariate=F)$psrf[,2],gelman.diag(m_tot_VCV,multivariate=F)$psrf[,2]) >= 1.1) #ok

#Save outputs
write.table(aa.prop,"Coefficients_proportions_MCMCglmm.txt",sep="\t")

