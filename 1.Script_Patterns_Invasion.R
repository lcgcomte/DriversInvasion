rm(list=ls(all=TRUE))
#/////////////////////////////////////////
#calculate life-history affinities
#//////////////////////////////////////////

#Occurrence database at the HUC8 scale (clean with names harmonized)
dat_fin= read.table("database_Fish.txt",h = T)  

#Trait database (clean with imputed values)
life = read.table("database_traits_imputed.txt",h=T) 

#ln-transformation traits
life$FECUN = log(life$FECUN)
life$AGEMAT = log(life$AGEMAT)

#Scale traits between 0 and 1
life$FECUN = (life$FECUN-min(life$FECUN,na.rm=T))/(max(life$FECUN,na.rm=T)-min(life$FECUN,na.rm=T)) 
life$AGEMAT = (life$AGEMAT-min(life$AGEMAT,na.rm=T))/(max(life$AGEMAT,na.rm=T)-min(life$AGEMAT,na.rm=T)) 
life$PC = (life$PC-min(life$PC,na.rm=T))/(max(life$PC,na.rm=T)-min(life$PC,na.rm=T)) 

#Creating endpoints life-history triangle
equilibrium = data.frame(FECUN = 0.5,PC = 1,AGEMAT = 1)
opportunistic = data.frame(FECUN = 0,PC = 0,AGEMAT =0)
periodic = data.frame(FECUN = 1,PC = 0,AGEMAT = 1)

#Calculating inverse distance between species and endpoints (=affinities)
require(cluster);

affinities = NULL
for(i in 1:nrow(life)){
dat = rbind(life[i,2:4],equilibrium)
affinities.equi <- 1 / daisy(dat,metric="euclidean")[1]

dat = rbind(life[i,2:4],opportunistic)
affinities.oppo <- 1 / daisy(dat,metric="euclidean")[1]

dat = rbind(life[i,2:4],periodic)
affinities.per <- 1 / daisy(dat,metric="euclidean")[1]

affinities = rbind(affinities,c(affinities.equi,affinities.oppo,affinities.per))
}

rownames(affinities) = life$Species_name
affinities_fin = (affinities/apply(affinities,1,sum))
colnames(affinities_fin) = c("equilibrium","opportunistic","periodic")
write.table(affinities_fin,"triangle_life-history.txt",sep="\t")

rm(list=ls(all=TRUE))
#/////////////////////////////////////////////////////////////////////
#Compute ln-ratios
#////////////////////////////////////////////////////////////////////
require(dataRetrieval);

affinities_fin = read.table("triangle_life-history.txt",h = T)
dat_fin = read.table("database_Fish.txt",h = T)  
dat_fin$HUC8 = zeroPad(dat_fin$HUC8,8)  #pad HUC8 with zeros

#File with characteristics of HUC8 [available at http://nhd.usgs.gov/wbd.html]
HUC = read.table("WBD_HUC8.txt",h=T)
HUC$HUC8 = zeroPad(HUC$HUC8,8)  #pad HUC8 with zeros

NL = levels(factor(dat_fin$HUC8))

#Compute richness life-history strategies for nonnative only
strategies_NAS = NULL
for(i in 1:length(NL)){
dat = dat_fin[dat_fin$HUC8 == NL[i] & dat_fin$ORIGIN != "Native",]
strategies_NAS = rbind(strategies_NAS,apply(affinities_fin[rownames(affinities_fin) %in% dat$Species_name,1:3],2,sum,na.rm=T))
}
rownames(strategies_NAS) = NL

#Compute richness life-history strategies for native only
strategies_NA = NULL
for(i in 1:length(NL)){
dat = dat_fin[dat_fin$HUC8 == NL[i] & dat_fin$ORIGIN == "Native",]
strategies_NA = rbind(strategies_NA,apply(affinities_fin[rownames(affinities_fin) %in% dat$Species_name,1:3],2,sum,na.rm=T))
}
rownames(strategies_NA) = NL

#Compute proportions life-history strategies
prop_NA = strategies_NA/apply(strategies_NA,1,sum)
prop_NAS = strategies_NAS/apply(strategies_NAS,1,sum)
prop_NA[is.na(prop_NA)] = 0
prop_NAS[is.na(prop_NAS)] = 0

#Compute species densities
Area = HUC$AreaSqKm[match(rownames(strategies_NAS),HUC$HUC8)] 
dens_NAS = apply(strategies_NAS,1,sum)/Area
dens_NA = apply(strategies_NA,1,sum)/Area

#ln-ratios (to NA)
lnratio_prop = log((prop_NAS+1)/(prop_NA+1))
lnratio_dens = log((dens_NAS+0.0001)/(dens_NA+0.0001))

#Save files
write.table(lnratio_prop,"LnRatio_proportion_life-history_NAS_NA.txt",sep="\t")
write.table(lnratio_dens,"LnRatio_densities_NAS_NA.txt",sep="\t")
write.table(dens_NA,"Densities_Natives.txt",sep="\t")

#/////////////////////////////////////////////////////////////////////
#Compute functional overlap
#/////////////////////////////////////////////////////////////////////
dist_NA = NULL
for(i in 1:length(NL)){

#Identify native species pool
dat.NA = dat_fin[dat_fin$HUC8 == NL[i] & dat_fin$ORIGIN == "Native",]  
affinities.NA = affinities_fin[rownames(affinities_fin) %in% dat.NA$Species_name,1:3]

#Identify nonnative species pool
dat.NAS = dat_fin[dat_fin$HUC8 == NL[i] & dat_fin$ORIGIN != "Native",] 
affinities.NAS = affinities_fin[rownames(affinities_fin) %in% dat.NAS$Species_name,1:3]

#Compute centroids
NA.cent = apply(affinities.NA,2,mean,na.rm=T)
NAS.cent = apply(affinities.NAS,2,mean,na.rm=T)

#Calculate Euclidean distance
dist_NA = rbind(dist_NA,dist(rbind(NA.cent,NAS.cent)))
}
rownames(dist_NA) = NL

#Take inverse distance (=overlap)
dist_NA = 1/dist_NA
dist_NA[is.na(dist_NA)] = 0   #set overlap to zero when only one group

#Save file
write.table(dist_NA,"Overlap_NAS_NA.txt",sep="\t")

rm(list=ls(all=TRUE))
#/////////////////////////////////////////////////////////////////////
#Assess regional trends (HUC2 hydrologic unit)
#/////////////////////////////////////////////////////////////////////
lnratio_dens = read.table("LnRatio_densities_NAS_NA.txt",h = T)
lnratio_prop = read.table("LnRatio_proportion_life-history_NAS_NA.txt",h = T)

HUC2 = sapply(strsplit(rownames(lnratio_prop),"*"),function(x) paste0(x[1],x[2],sep=""))

#///////////////////////////
#ln-ratio densities
#//////////////////////////
#Linear model without a separate intercept for each hydrologic unit
m.dens = lm(lnratio_dens$x ~ HUC2-1)
anova(m.dens)
summary(m.dens)

#///////////////////////////
#ln-ratio proportions
#//////////////////////////
#equilibrium
m.eq = lm(lnratio_prop$equilibrium ~ HUC2-1)
anova(m.eq)
summary(m.eq)

#opportunistic
m.op = lm(lnratio_prop$opportunistic ~ HUC2-1)
anova(m.op)
summary(m.op)

#periodic
m.pe = lm(lnratio_prop$periodic ~ HUC2-1)
anova(m.pe)
summary(m.pe)

rm(list=ls(all=TRUE))
#////////////////////////////////////////////////////////////////////
#Null Model random colonizations 
#///////////////////////////////////////////////////////////////////
require(dataRetrieval);

affinities_fin = read.table("triangle_life-history.txt",h = T)
dat_fin = read.table("database_Fish.txt",h = T)  
dat_fin$HUC8 = zeroPad(dat_fin$HUC8,8)  #pad HUC8 with zeros

#File with characteristics of HUC8 [available at http://nhd.usgs.gov/wbd.html]
HUC = read.table("WBD_HUC8.txt",h=T)
HUC$HUC8 = zeroPad(HUC$HUC8,8)  #pad HUC8 with zeros

NL = levels(factor(dat_fin$HUC8))

dat_NAS = dat_fin[which(! dat_fin$ORIGIN == "Native"),]
dat_NA = dat_fin[which(dat_fin$ORIGIN == "Native"),]

#Species prevalence entire dataset
SP = table(factor(dat_NAS$Species_name))

dens_NAS_null = prop_NAS_null = NULL
for(j in 1:999){
##reassign species at random
dat_null = NULL
for(i in 1:length(SP)){
nat = dat_NA$HUC8[dat_NA$New_name == names(SP[i])]   #remove basins where species is native
if(length(nat) > 0){
NL_NAS = NL[-which(NL %in% nat)]
} else{
NL_NAS = NL
}
dat_null = rbind(dat_null, cbind(names(SP[i]),sample(NL_NAS,SP[i])))
}
dat_null = data.frame(dat_null)
names(dat_null) = c("Species_name","HUC8")

#SR per huc8 for NAS only
strategies_NAS = NULL
for(i in 1:length(NL)){
dat = dat_null[dat_null$HUC8 == NL[i] ,]
strategies_NAS = rbind(strategies_NAS,apply(affinities_fin[rownames(affinities_fin) %in% dat$Species_name,1:3],2,sum,na.rm=T))
}
rownames(strategies_NAS) = NL               

#Proportions
prop_NAS = strategies_NAS/apply(strategies_NAS,1,sum)
prop_NAS[is.na(prop_NAS)] = 0

#Densities
Area = HUC$AreaSqKm[match(rownames(strategies_NAS),HUC$HUC8)] 
dens_NAS = apply(strategies_NAS,1,sum)/Area

#SR per huc8 for NA only
if(j == 1){ #will always be the same so better to calculate only once
strategies_NA = NULL
for(i in 1:length(NL)){
dat = dat_fin[dat_fin$HUC8 == NL[i] & dat_fin$ORIGIN == "Native",]
strategies_NA = rbind(strategies_NA,apply(affinities_fin[rownames(affinities_fin) %in% dat$Species_name,1:3],2,sum,na.rm=T))
}
rownames(strategies_NA) = NL

dens_NA = apply(strategies_NA,1,sum)/Area

prop_NA = strategies_NA/apply(strategies_NA,1,sum)
prop_NA[is.na(prop_NA)] = 0
}

#ln-ratios (to NA)
lnratio_prop = log((prop_NAS+1)/(prop_NA+1))
lnratio_dens = log((dens_NAS+0.0001)/(dens_NA+0.0001))

#save results
prop_NAS_null = cbind(prop_NAS_null,lnratio_prop)
dens_NAS_null = cbind(dens_NAS_null,lnratio_dens)

print(j)
}

write.table(dens_NAS_null,"lnRatio_densities_NAS_null_model.txt",sep="\t")
write.table(prop_NAS_null,"lnRatio_proportions_NAS_null_model.txt",sep="\t")

rm(list=ls(all=TRUE))
#/////////////////////////////////
#Comparaison with observed values
#/////////////////////////////////
#observed
lnratio_prop=read.table("LnRatio_proportion_life-history_NAS_NA.txt",h=T)
lnratio_dens=read.table("LnRatio_densities_NAS_NA.txt",h=T)

#predicted
lnratio_dens_null = read.table("lnRatio_densities_NAS_null_model.txt",h=T)
lnratio_prop_null = read.table("lnRatio_proportions_NAS_null_model.txt",h=T)

Pval = Pval.eq = Pval.op = Pval.pe = NULL
for(i in 1:nrow(lnratio_prop_null)){

#proportion of simulations inferior or superior to observed value
tail1 = length(which(lnratio_dens_null[i,] <= lnratio_dens[i,1]))/999
tail2 = length(which(lnratio_dens_null[i,] >= lnratio_dens[i,1]))/999
Pval = c(Pval,ifelse(tail1 < 0.05 , "lower",ifelse(tail2 < 0.05,"higher","non-significant")))

ratio_prop = lnratio_prop_null[i,seq(1,ncol(lnratio_prop_null),3)]
tail1.eq = length(which(ratio_prop <= lnratio_prop[i,1]))
tail2.eq = length(which(ratio_prop >= lnratio_prop[i,1]))
Pval.eq = c(Pval.eq,ifelse(tail1.eq/999 < 0.05,"lower",ifelse(tail2.eq/999 < 0.05,"higher","non-significant")))

ratio_prop = lnratio_prop_null[i,seq(2,ncol(lnratio_prop_null),3)]
tail1.op = length(which(ratio_prop <= lnratio_prop[i,2]))
tail2.op = length(which(ratio_prop >= lnratio_prop[i,2]))
Pval.op = c(Pval.op,ifelse(tail1.op/999 < 0.05,"lower",ifelse(tail2.op/999 < 0.05,"higher","non-significant")))

ratio_prop = lnratio_prop_null[i,seq(3,ncol(lnratio_prop_null),3)]
tail1.pe = length(which(ratio_prop <= lnratio_prop[i,3]))
tail2.pe = length(which(ratio_prop >= lnratio_prop[i,3]))
Pval.pe = c(Pval.pe,ifelse(tail1.pe/999 < 0.05,"lower",ifelse(tail2.pe/999 < 0.05,"higher","non-significant")))

print(i)
}

res = data.frame(Pval,Pval.eq,Pval.op,Pval.pe)
rownames(res) = rownames(lnratio_prop_null)

write.table(res,"Results_null_model.txt",sep="\t")
