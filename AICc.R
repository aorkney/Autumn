# Compute AICc values for different GLS models 
# among variables in the Barents Sea Septembers dataset
# You will need to download, from the Github directory, 
# the dataset 'All_variables_fall_21_08_2021.csv' 
# which contains the gridded data for analysis 


years<-c(2002:2019)

library(nlme) # package for fitting GLS models
# if you do not have this, run 'install.packages('nlme')'

setwd('D:/Documents/OCCCI2/Data_store')
# You will need to adjust your work directory accordingly 

dataset<-read.csv('All_variables_fall_21_08_2021.csv')
# load the re-gridded data for SST, [chl-a], winds, currents, salinity

deg2rad <- function(deg) {(deg * pi) / (180)} # function to convert degrees to radians 

current_ts<-list()
for(j in years){
	current_ts[[j]]<-weighted.mean(x=dataset$velocity[which(dataset$lat<75 & dataset$lon<20 & dataset$years==j)],na.rm=T,
	w= cos(deg2rad(dataset$lat[which(dataset$lat<75 & dataset$lon<20 & dataset$years==j)])) )
}
current_ts<-unlist(current_ts)
# Time series of current velocity weighted by surface area in the Atlantic inflow region


biomass_ts<-list()
for(j in years){
	biomass_ts[[j]]<-weighted.mean(x=dataset$chla[which(dataset$years==j)],na.rm=T,
	w= cos(deg2rad(dataset$lat[which(dataset$years==j)])) )
}
biomass_ts<-unlist(biomass_ts)
# Time series of [chl-a] weighted by surface area 

sst_ts<-list()
for(j in years){
	sst_ts[[j]]<-weighted.mean(x=dataset$sst[which(dataset$years==j)],na.rm=T,
	w= cos(deg2rad(dataset$lat[which(dataset$years==j)])) )
}
sst_ts<-unlist(sst_ts)
# Time series of SST weighted by surface area 

wind_ts<-list()
for(j in years){
	wind_ts[[j]]<-weighted.mean(x=dataset$storminess[which(dataset$years==j)][complete.cases(dataset$storminess[which(dataset$years==j)])],
	w= cos(deg2rad(dataset$lat[which(dataset$years==j)][complete.cases(dataset$storminess[which(dataset$years==j)])])) )
}
wind_ts<-unlist(wind_ts)
# Time series of windye vent frequency, weighted by surface area 

bio_nV<-(biomass_ts/sd(biomass_ts))-mean(biomass_ts/sd(biomass_ts))
sst_nV<-(sst_ts/sd(sst_ts))-mean(sst_ts/sd(sst_ts))
wind_nV<-(wind_ts/sd(wind_ts))-mean(wind_ts/sd(wind_ts))
current_nV<-(current_ts/sd(current_ts))-mean(current_ts/sd(current_ts))
# variance normalised time series 


df1<-as.data.frame(cbind(bio_nV,current_nV,years)) # example
# if you wish to explore other combinations of variables, 
# change the arguments being placed in the dataframe

colnames(df1)<-c('V1','V2','V3')
mod1<-gls(V1 ~ V2,data=df1) # standard model 
mod2<-gls( V1 ~ V2  ,correlation=corAR1(form=~V3),data= df1) # temporal autocorrelation
anova(mod1,mod2)$p[2] # which model is preferred? (p<0.05 indicates autocorrelated model is favoured)



library(AICcmodavg) # for computing corrected AIC statistic
# if you do not have this package run 'install.packages('AICcmodavg')'

AICc(mod1) # what is the AIC of the model?
# Change the model being tested by changing the argument to AICc 

# review other model properties by typing 'summary(mod1)'


