# This script will produce plots of the linear relationships between
# [chl-a] inferred from the OC-CCI record: 
# https://www.oceancolour.org/thredds/catalog-cci.html
# using the OC6_MERIS algorithm of O'reilly & Werdell, 2019 
# and average current velocity in the Atlantic inflow region, inferred from 
# ARCTIC REANALYSIS PHYS 002 003 available on https://marine.copernicus.eu/
# the inter-annual differences of these time series;
# and between the inter-annual differences of [chl-a] and
# Windy events, inferred from WIND GLO WIND L3 REP OBSERVATIONS 012 005
# available on https://marine.copernicus.eu/
# Also displayed are variance-normalised time series of 
# [chl-a], SST, inflowing current velocity and windy event frequency, 
# at a basin scale. 
# Contours of Atlantic watermass distribution are projected on maps
# to foster clearer interpretation.
# SST and salinity data are vailable from
# ARCTIC REANALYSIS PHYS 002 003 available on https://marine.copernicus.eu/
# This analysis was coded in R v 3.63 by Andrew Orkney
# Storminess in the average number of days with average winds above 10 metres per second that
# are observed in a monthly interval 

# You will need to download, from the Github directory, 
# the bathymetric dataset 'marmap_coord_0;70;80;85_res_5.csv'
# which was constructed from the ETOPO1 NOAA bathymetric database
# and the dataset 'All_variables_fall_21_08_2021.csv' 
# which contains the gridded data for analysis 

# Load Barents bathymetry and re-project in laea
setwd("D:/Documents/cruise_4")
# You will need to adjust the working directory accordingly

library(marmap)# For projecting maps 
# You will need to run 'install.packages('marmap') if you do not have this package

Barents<-read.csv("marmap_coord_0;70;80;85_res_5.csv")# load data
Barents<-Barents[-which(Barents$V1>60),]# exclude data in the far east
Barents<-as.bathy(Barents)# convert to bathymetric object
r1<-as.raster(Barents)# rasterise
projection<-"+proj=laea +lat_0=70, +lon_0=30" # define a projection
library(raster) # For raster projections
# If you do not have this package, install it as before
r2<-projectRaster(r1,crs=projection) # re-project the data
as.bathy(r2) ->bath2 # convert back to a bathymetric object

# trouble shooting: detach("package:marmap",unload=T)
# you may need to detach and reload the package 'marmap'


fortify.bathy(bath2)->bath3
bath3$z[which(bath3$z>=0)]<-0
bath4<-bath3
bath4$z[which(bath4$z<0)]<-NA
# This has produced a land mask 

label_pos<-cbind(c(0,30,60,66,68),c(69.5,69.5,69.5,75,80))
d<-data.frame(label_pos)
colnames(d)<-c('lon','lat')
coordinates(d)<-c("lon","lat")
library(rgdal)
proj4string(d)<-CRS("+init=epsg:4326")
CRS.new<-CRS("+proj=laea +lat_0=70, +lon_0=30")
d.new<-spTransform(d,CRS.new)
# This is just to produce labels in the final plot for lat-lon coords 


# Now we need to run analyses

setwd('D:/Documents/OCCCI2/Data_store')
# You will need to adjust your work directory accordingly 


dataset<-read.csv('All_variables_fall_21_08_2021.csv')
# load the re-gridded data for SST, [chl-a], winds, currents, salinity

dataset<-cbind(dataset,NA)
# Add a new blank column 

colnames(dataset)[15]<-'Atlantic'
dataset$Atlantic<-0
dataset$Atlantic[which(dataset$sst>=3 & dataset$salinity >=34.8)]<-1
# Define a column for defining the water-mass based on hydrography
# we may or may not need this in downstream analysis 

years<-c(2002:2019)
lat_steps<-seq(70,85,1)
lon_steps<-seq(0,60,1)
# Defining temporal and spatial steps 

# You need to make a time series of the average current speed evolution south of 75N and west of 20E

# You need to make sure that the current speed dataset is 
# weighted according to its latitude 

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



trends<-cbind(expand.grid(lat_steps,lon_steps),NA,NA)
trends<-as.data.frame(trends)
trends<-cbind(trends,rep(0,dim(trends)[1]))
trends<-cbind(trends,rep(0,dim(trends)[1]))
colnames(trends)<-c('lat','lon','trends','odds','auto','Atlantic')

library(nlme) # this is a package for fitting regressions
# if you do not have it, install it as before

# we are going to cycle over unique lat-lon combinations
# and compute trends now 


for(a in 1:(length(lat_steps)-1)){
	for(b in 1:(length(lon_steps)-1)){
		waters<-dataset[which(dataset$lat==lat_steps[a] & dataset$lon==lon_steps[b]),]$Atlantic
			if(2>1){
			chosen_row<-which(trends$lat==lat_steps[a] & trends$lon==lon_steps[b])
			trends[chosen_row,6]<-mean(waters,na.rm=T)
			}
	}
}
# Atlantic watermass distribution 

detach("package:marmap",unload=T)
library(marmap)
contour<-as.bathy(trends[,c(2,1,6)])
contour_raster<-as.raster(contour)
contour_2<-projectRaster(contour_raster,crs=projection)
as.bathy(contour_2)->ready_contour
fortify.bathy(ready_contour)->contour_3

# Contours of Atlantic watermass distribution 


library(nlme)# this is a package for fitting regressions
# if you do not have it, install it as before

# we are going to cycle over unique lat-lon combinations
# and compute trends now 


for(a in 1:(length(lat_steps)-1)){ # for each latitude step 
	for(b in 1:(length(lon_steps)-1)){ # for each longitude step 
		biomass<-dataset[which(dataset$lat==lat_steps[a] & dataset$lon==lon_steps[b]),]$chla # fetch [chl-a] data 
			if( length(which(is.na(current_ts)==F))>3 & sd(current_ts,na.rm=T)>0 & length(which(is.na(biomass)==F))>3 & sd(biomass,na.rm=T)>0){ # check data good 
			chosen_row<-which(trends$lat==lat_steps[a] & trends$lon==lon_steps[b])
			df1<-as.data.frame(cbind(biomass,current_ts,years))
			df1<-df1[complete.cases(df1),]
			mod1<-gls( biomass ~ current_ts ,data=df1) # regression of biomass and current speed 
			tryCatch( mod2<-gls( biomass ~ current_ts  ,correlation=corAR1(form=~years),data= df1) ,error=function(e) print('test') )
			# fit an autocorrelated model as well 
			pval<-1
			if(exists('mod2')){
					pval<-anova(mod1,mod2)$p[2] # which model is preferred? 
			} 

			if(pval <0.05){
					trends[chosen_row,3]<-summary(mod2)$t[2,1]
					trends[chosen_row,4]<-summary(mod2)$t[2,4]
					trends[chosen_row,5]<-2
					rm('mod2')	
				
			} else {
				trends[chosen_row,3]<-summary(mod1)$t[2,1]
				trends[chosen_row,4]<-summary(mod1)$t[2,4]
				trends[chosen_row,5]<-1
				if(exists('mod2')){
					rm(mod2)
				}
				
			}
		
			
			}
	}
	print(a) # print latitude step (1:15) to update on progress
}

library(pracma)# This is a package we need for the regridding of the results
# if you do not have it, install it as before 


xp <- linspace(min(lon_steps),max(lon_steps),150)
yp <- linspace(min(lat_steps),max(lat_steps),150)
# re-grid 

mgrid<-meshgrid(yp,xp)
Z <- 0* mgrid$X

for(i in 1:(dim(trends)[1]-1) ){
	Z[which(mgrid$Y >= trends$lon[i] & mgrid$Y <= (trends$lon[i]+1)
	& mgrid$X >= trends$lat[i] & mgrid$X <= (trends$lat[i]+1))]<-trends$trends[i]

}

newdat<-as.data.frame(cbind(expand.grid(xp,yp),as.numeric(Z)))
colnames(newdat)<-c('lon','lat','trends')
# re-gridding is complete 

detach("package:marmap",unload=T)
library(marmap)

# Now we can re-project the results 
Chlor<-as.bathy(newdat[,c(1,2,3)])
Chlor_raster<-as.raster(Chlor)
Chlor_2<-projectRaster(Chlor_raster,crs=projection)
as.bathy(Chlor_2)->ready_Chlor
fortify.bathy(ready_Chlor)->Chlor_3
# 'Chlor_3' is a reprojected summary of the relationship of [chl-a] to inflow current speed 

probs<-as.bathy(trends[,c(2,1,4)])
probs_raster<-as.raster(probs)
probs_2<-projectRaster(probs_raster,crs=projection)
as.bathy(probs_2)->ready_probs
fortify.bathy(ready_probs)->probs_3
# 'probs_3' is a summary of the p-valus associated with the preferred model fits 

stip_pos<-trends[which(trends$odds<=0.05),c(2,1)]+0.5
stip_d<-data.frame(stip_pos)
colnames(stip_d)<-c('lon','lat')
coordinates(stip_d)<-c("lon","lat")
proj4string(stip_d)<-CRS("+init=epsg:4326")
CRS.new<-CRS("+proj=laea +lat_0=70, +lon_0=30")
stip_d.new<-spTransform(stip_d,CRS.new)
point_shape<-trends[which(trends$odds<=0.05),5]
# This is a vector specifying whether standard or autocorrelated fits are preferred 

hard<-max(abs(trends$trends)[which(trends$odds<0.05)],na.rm=T)
limits<-c(-hard,hard)
Chlor_3$z[which(Chlor_3$z< -hard)]<- -hard
Chlor_3$z[which(Chlor_3$z> hard)]<- hard
# define scale limits

library(ggplot2)# package for producing aesthetic plots
# if you do not have this, install it as before 

B_v_Ci<-
ggplot(data=Chlor_3, aes(x=x,y=y))+
coord_equal()+
geom_raster(aes(x=x,y=y,fill=z),interpolate=F)+
scale_fill_gradient2(high='red',low='blue',mid='white',limits=limits,na.value='white')+
geom_contour(data=contour_3,aes(x=x,y=y,z=z),breaks=.05,
col='black')+
geom_contour(data=contour_3,aes(x=x,y=y,z=z),breaks=.5,
col='black',size=0.75)+
geom_contour(data=contour_3,aes(x=x,y=y,z=z),breaks=.95,
col='black',size=1.5)+
geom_raster(data = dplyr::filter(bath4, !is.na(z)), 
            aes(x = x, y = y), fill = "gray40")+
geom_point(data=as.data.frame(coordinates(stip_d.new)),
aes(x=lon,y=lat),alpha=.5,size=1.2,shape=point_shape)+
geom_text(data=as.data.frame(coordinates(d.new)),
aes(x=lon,y=lat), label=c(expression(paste('0'^o,'E')),expression(paste('30'^o,'E')),expression(paste('60'^o,'E')),expression(paste('75'^o,'N')),expression(paste('80'^o,'N'))),
angle=c(-30,0,30,30,30))+
labs(fill=expression(paste(Delta, italic('B C'[i]),''^-1)) )+
theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="bottom",
legend.title=element_text(size=16),
legend.text=element_text(size=10),
legend.key.width = unit(3, "line"),
legend.key.height = unit(2, "line"),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
# type 'B_v_Ci' into the terminal to view the plot
# significant trends at a 0.05 alpha level in [chl-a] against inflow current speed are marked with a circle (standard fit) or triangle (autocorrelated)

# We are going to now perform the same exercise for the remaining variables
# I will desist comments for now 

trends<-cbind(expand.grid(lat_steps,lon_steps),NA,NA)
trends<-as.data.frame(trends)
trends<-cbind(trends,rep(0,dim(trends)[1]))
trends<-cbind(trends,rep(0,dim(trends)[1]))
colnames(trends)<-c('lat','lon','trends','odds','auto','Atlantic')


for(a in 1:(length(lat_steps)-1)){
	for(b in 1:(length(lon_steps)-1)){
		waters<-dataset[which(dataset$lat==lat_steps[a] & dataset$lon==lon_steps[b]),]$Atlantic
			if(2>1){
			chosen_row<-which(trends$lat==lat_steps[a] & trends$lon==lon_steps[b])
			trends[chosen_row,6]<-mean(waters,na.rm=T)
			}
	}
}

contour<-as.bathy(trends[,c(2,1,6)])
contour_raster<-as.raster(contour)
contour_2<-projectRaster(contour_raster,crs=projection)
as.bathy(contour_2)->ready_contour
fortify.bathy(ready_contour)->contour_3



for(a in 1:(length(lat_steps)-1)){
	for(b in 1:(length(lon_steps)-1)){
		biomass<-dataset[which(dataset$lat==lat_steps[a] & dataset$lon==lon_steps[b]),]$chla
			if( length(which(is.na(diff(current_ts))==F))>3 & sd(diff(current_ts),na.rm=T)>0 & length(which(is.na(diff(biomass))==F))>3 & sd(diff(biomass),na.rm=T)>0){
			chosen_row<-which(trends$lat==lat_steps[a] & trends$lon==lon_steps[b])
			df1<-as.data.frame(cbind(diff(biomass),diff(current_ts),years[-1]))
			df1<-df1[complete.cases(df1),]
			mod1<-gls( V1 ~ V2 ,data=df1)
			tryCatch( mod2<-gls( V1 ~ V2  ,correlation=corAR1(form=~V3),data= df1) ,error=function(e) print('test') )
			pval<-1
			if(exists('mod2')){
					pval<-anova(mod1,mod2)$p[2]
			} 

			if(pval <0.05){
					trends[chosen_row,3]<-summary(mod2)$t[2,1]
					trends[chosen_row,4]<-summary(mod2)$t[2,4]
					trends[chosen_row,5]<-2
					rm('mod2')	
				
			} else {
				trends[chosen_row,3]<-summary(mod1)$t[2,1]
				trends[chosen_row,4]<-summary(mod1)$t[2,4]
				trends[chosen_row,5]<-1
				if(exists('mod2')){
					rm(mod2)
				}
				
			}
		
			
			}
	}
	print(a)
}



xp <- linspace(min(lon_steps),max(lon_steps),150)
yp <- linspace(min(lat_steps),max(lat_steps),150)

mgrid<-meshgrid(yp,xp)
Z <- 0* mgrid$X

for(i in 1:(dim(trends)[1]-1) ){
	Z[which(mgrid$Y >= trends$lon[i] & mgrid$Y <= (trends$lon[i]+1)
	& mgrid$X >= trends$lat[i] & mgrid$X <= (trends$lat[i]+1))]<-trends$trends[i]

}

newdat<-as.data.frame(cbind(expand.grid(xp,yp),as.numeric(Z)))
colnames(newdat)<-c('lon','lat','trends')


# Now we can re-project the results 
Chlor<-as.bathy(newdat[,c(1,2,3)])
Chlor_raster<-as.raster(Chlor)
Chlor_2<-projectRaster(Chlor_raster,crs=projection)
as.bathy(Chlor_2)->ready_Chlor
fortify.bathy(ready_Chlor)->Chlor_3

probs<-as.bathy(trends[,c(2,1,4)])
probs_raster<-as.raster(probs)
probs_2<-projectRaster(probs_raster,crs=projection)
as.bathy(probs_2)->ready_probs
fortify.bathy(ready_probs)->probs_3

stip_pos<-trends[which(trends$odds<=0.05),c(2,1)]+0.5
stip_d<-data.frame(stip_pos)
colnames(stip_d)<-c('lon','lat')
coordinates(stip_d)<-c("lon","lat")
proj4string(stip_d)<-CRS("+init=epsg:4326")
CRS.new<-CRS("+proj=laea +lat_0=70, +lon_0=30")
stip_d.new<-spTransform(stip_d,CRS.new)
point_shape<-trends[which(trends$odds<=0.05),5]

hard<-max(abs(trends$trends)[which(trends$odds<0.05)],na.rm=T)
limits<-c(-hard,hard)
Chlor_3$z[which(Chlor_3$z< -hard)]<- -hard
Chlor_3$z[which(Chlor_3$z> hard)]<- hard

dB_v_dCi<-
ggplot(data=Chlor_3, aes(x=x,y=y))+
coord_equal()+
geom_raster(aes(x=x,y=y,fill=z),interpolate=F)+
scale_fill_gradient2(high='red',low='blue',mid='white',limits=limits,na.value='white')+
geom_contour(data=contour_3,aes(x=x,y=y,z=z),breaks=.05,
col='black')+
geom_contour(data=contour_3,aes(x=x,y=y,z=z),breaks=.5,
col='black',size=0.75)+
geom_contour(data=contour_3,aes(x=x,y=y,z=z),breaks=.95,
col='black',size=1.5)+
geom_raster(data = dplyr::filter(bath4, !is.na(z)), 
            aes(x = x, y = y), fill = "gray40")+
geom_point(data=as.data.frame(coordinates(stip_d.new)),
aes(x=lon,y=lat),alpha=.5,size=1.2,shape=point_shape)+
geom_text(data=as.data.frame(coordinates(d.new)),
aes(x=lon,y=lat), label=c(expression(paste('0'^o,'E')),expression(paste('30'^o,'E')),expression(paste('60'^o,'E')),expression(paste('75'^o,'N')),expression(paste('80'^o,'N'))),
angle=c(-30,0,30,30,30))+
labs(fill=expression( paste(Delta, italic(dot(B)),' ',italic(dot(C)[i]))^-1 ))+
theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="bottom",
legend.title=element_text(size=16),
legend.text=element_text(size=10),
legend.key.width = unit(3, "line"),
legend.key.height = unit(2, "line"),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())


trends<-cbind(expand.grid(lat_steps,lon_steps),NA,NA)
trends<-as.data.frame(trends)
trends<-cbind(trends,rep(0,dim(trends)[1]))
trends<-cbind(trends,rep(0,dim(trends)[1]))
colnames(trends)<-c('lat','lon','trends','odds','auto','Atlantic')


for(a in 1:(length(lat_steps)-1)){
	for(b in 1:(length(lon_steps)-1)){
		waters<-dataset[which(dataset$lat==lat_steps[a] & dataset$lon==lon_steps[b]),]$Atlantic
			if(2>1){
			chosen_row<-which(trends$lat==lat_steps[a] & trends$lon==lon_steps[b])
			trends[chosen_row,6]<-mean(waters,na.rm=T)
			}
	}
}

contour<-as.bathy(trends[,c(2,1,6)])
contour_raster<-as.raster(contour)
contour_2<-projectRaster(contour_raster,crs=projection)
as.bathy(contour_2)->ready_contour
fortify.bathy(ready_contour)->contour_3



for(a in 1:(length(lat_steps)-1)){
	for(b in 1:(length(lon_steps)-1)){
		biomass<-dataset[which(dataset$lat==lat_steps[a] & dataset$lon==lon_steps[b]),]$chla
		wind<-dataset[which(dataset$lat==lat_steps[a] & dataset$lon==lon_steps[b]),]$storminess
			if( length(which(is.na(diff(wind))==F))>3 & sd(diff(wind),na.rm=T)>0 & length(which(is.na(diff(biomass))==F))>3 & sd(diff(biomass),na.rm=T)>0){
			chosen_row<-which(trends$lat==lat_steps[a] & trends$lon==lon_steps[b])
			df1<-as.data.frame(cbind(diff(biomass),diff(wind),years[-1]))
			df1<-df1[complete.cases(df1),]
			if(dim(df1)[1] >2){
			mod1<-gls( V1 ~ V2 ,data=df1)
			tryCatch( mod2<-gls( V1 ~ V2  ,correlation=corAR1(form=~V3),data= df1) ,error=function(e) print('test') )
			pval<-1
			if(exists('mod2')){
					pval<-anova(mod1,mod2)$p[2]
			} 

			if(pval <0.05){
					trends[chosen_row,3]<-summary(mod2)$t[2,1]
					trends[chosen_row,4]<-summary(mod2)$t[2,4]
					trends[chosen_row,5]<-2
					rm('mod2')	
				
			} else {
				trends[chosen_row,3]<-summary(mod1)$t[2,1]
				trends[chosen_row,4]<-summary(mod1)$t[2,4]
				trends[chosen_row,5]<-1
				if(exists('mod2')){
					rm(mod2)
				}
				
			}
		
			
			}
	}
}
	print(a)
}



xp <- linspace(min(lon_steps),max(lon_steps),150)
yp <- linspace(min(lat_steps),max(lat_steps),150)

mgrid<-meshgrid(yp,xp)
Z <- 0* mgrid$X

for(i in 1:(dim(trends)[1]-1) ){
	Z[which(mgrid$Y >= trends$lon[i] & mgrid$Y <= (trends$lon[i]+1)
	& mgrid$X >= trends$lat[i] & mgrid$X <= (trends$lat[i]+1))]<-trends$trends[i]

}

newdat<-as.data.frame(cbind(expand.grid(xp,yp),as.numeric(Z)))
colnames(newdat)<-c('lon','lat','trends')


# Now we can re-project the results 
Chlor<-as.bathy(newdat[,c(1,2,3)])
Chlor_raster<-as.raster(Chlor)
Chlor_2<-projectRaster(Chlor_raster,crs=projection)
as.bathy(Chlor_2)->ready_Chlor
fortify.bathy(ready_Chlor)->Chlor_3

probs<-as.bathy(trends[,c(2,1,4)])
probs_raster<-as.raster(probs)
probs_2<-projectRaster(probs_raster,crs=projection)
as.bathy(probs_2)->ready_probs
fortify.bathy(ready_probs)->probs_3

stip_pos<-trends[which(trends$odds<=0.05),c(2,1)]+0.5
stip_d<-data.frame(stip_pos)
colnames(stip_d)<-c('lon','lat')
coordinates(stip_d)<-c("lon","lat")
proj4string(stip_d)<-CRS("+init=epsg:4326")
CRS.new<-CRS("+proj=laea +lat_0=70, +lon_0=30")
stip_d.new<-spTransform(stip_d,CRS.new)
point_shape<-trends[which(trends$odds<=0.05),5]

hard<-max(abs(trends$trends)[which(trends$odds<0.05)],na.rm=T)
limits<-c(-hard,hard)
Chlor_3$z[which(Chlor_3$z< -hard)]<- -hard
Chlor_3$z[which(Chlor_3$z> hard)]<- hard

dB_v_dW<-
ggplot(data=Chlor_3, aes(x=x,y=y))+
coord_equal()+
geom_raster(aes(x=x,y=y,fill=z),interpolate=F)+
scale_fill_gradient2(high='red',low='blue',mid='white',limits=limits,na.value='white')+
geom_contour(data=contour_3,aes(x=x,y=y,z=z),breaks=.05,
col='black')+
geom_contour(data=contour_3,aes(x=x,y=y,z=z),breaks=.5,
col='black',size=0.75)+
geom_contour(data=contour_3,aes(x=x,y=y,z=z),breaks=.95,
col='black',size=1.5)+
geom_raster(data = dplyr::filter(bath4, !is.na(z)), 
            aes(x = x, y = y), fill = "gray40")+
geom_point(data=as.data.frame(coordinates(stip_d.new)),
aes(x=lon,y=lat),alpha=.5,size=1.2,shape=point_shape)+
geom_text(data=as.data.frame(coordinates(d.new)),
aes(x=lon,y=lat), label=c(expression(paste('0'^o,'E')),expression(paste('30'^o,'E')),expression(paste('60'^o,'E')),expression(paste('75'^o,'N')),expression(paste('80'^o,'N'))),
angle=c(-30,0,30,30,30))+
labs(fill=expression( paste(Delta, italic(dot(B)),' ',italic(dot(W)))^-1 ))+
theme(axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.position="bottom",
legend.title=element_text(size=16),
legend.text=element_text(size=10),
legend.key.width = unit(3, "line"),
legend.key.height = unit(2, "line"),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())

# Now the little time series subplot


bio_nV<-(biomass_ts/sd(biomass_ts))-mean(biomass_ts/sd(biomass_ts))
sst_nV<-(sst_ts/sd(sst_ts))-mean(sst_ts/sd(sst_ts))
wind_nV<-(wind_ts/sd(wind_ts))-mean(wind_ts/sd(wind_ts))
current_nV<-(current_ts/sd(current_ts))-mean(current_ts/sd(current_ts))
# variance normalised time series 

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# line colours 

linedata<-as.data.frame(cbind(years,bio_nV,sst_nV,current_nV,wind_nV))
# aggregate data into object

library(reshape2) # for re-organising data
# if you do not have this package then run 'install.packages('reshape2')'

long_line_data<-melt(linedata,id='years') # reorganised 

linetype=c(rep(1,18),rep(2,18),rep(4,18),rep(5,18))  # line style 

library(ggplot2) # plotting package (should be loaded already)

mp12<-
ggplot(data=long_line_data)+
geom_line(aes(x=years,y=value,colour=variable),size=0.6)+
geom_line(aes(x=years,y=value,colour=variable,linetype=variable),size=3/2 )+
	scale_colour_manual("",values=c(cbbPalette[4],cbbPalette[2],cbbPalette[1],cbbPalette[3],cbbPalette[6]), labels=c(expression(italic('B')),expression(italic('SST')),expression(italic('C'[i])),expression(italic('W')), expression(italic('Co')) )) +
  scale_linetype_manual("", values=c(1,2,4,5,6),labels=c(expression(italic('B')),expression(italic('SST')),expression(italic('C'[i])),expression(italic('W')), expression(italic('Co')) )) +
theme(axis.line=element_line(size=3/2),
      axis.text.x=element_text(size=14,colour='black',angle=-90,vjust=0.5,hjust=1),
      axis.text.y=element_text(size=14,colour='black'),
      axis.ticks=element_line(size=2),
      axis.title.x=element_text(size=10,colour='black'),
      axis.title.y=element_blank(),
      legend.position="bottom",
legend.title=element_text(size=10),
legend.text=element_text(size=10,colour=c('black'),face='bold'),
legend.key.width = unit(3, "line"),
legend.key.height = unit(2, "line"),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())+
  scale_x_continuous("years", labels = as.character(years), breaks = years)

library(ggpubr) # combine subplots
# if you do not have this package, install it as before 

dev.new(width=10,height=10,units='inches')
ggarrange(B_v_Ci,dB_v_dCi,dB_v_dW,mp12,ncol=2,nrow=2,labels=c('a','b','c','d'),label.x=0.175,label.y=1,font.label=list(size=24) )

#ggsave(filename='AtW_relationships_09_08_2021.pdf',dpi=300)
# uncomment this line if you wish to save the plot 
