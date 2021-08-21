# This script will produce plots of the average trends in
# [chl-a] inferred from the OC-CCI record: 
# https://www.oceancolour.org/thredds/catalog-cci.html
# using the OC6_MERIS algorithm of O'reilley & Werdell, 2019 
# Sea Surface Temperature, inferred from 
# ARCTIC REANALYSIS PHYS 002 003 available on https://marine.copernicus.eu/
# Windy events, inferred from WIND GLO WIND L3 REP OBSERVATIONS 012 005
# available on https://marine.copernicus.eu/
# and changes in current velocities, inferred from 
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

library(marmap) # For projecting maps 
# You will need to run 'install.packages('marmap') if you do not have this package

Barents<-read.csv("marmap_coord_0;70;80;85_res_5.csv") # load data
Barents<-Barents[-which(Barents$V1>60),] # exclude data in the far east
Barents<-as.bathy(Barents) # convert to bathymetric object
r1<-as.raster(Barents) # rasterise 
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

colnames(dataset)[10]<-'Atlantic'
dataset$Atlantic<-0
dataset$Atlantic[which(dataset$sst>=3 & dataset$salinity >=34.8)]<-1
# Define a column for defining the water-mass based on hydrography
# we may or may not need this in downstream analysis 

years<-c(2002:2019)
lat_steps<-seq(70,85,1)
lon_steps<-seq(0,60,1)
# Defining temporal and spatial steps 

trends<-cbind(expand.grid(lat_steps,lon_steps),NA,NA)
trends<-as.data.frame(trends)
trends<-cbind(trends,rep(0,dim(trends)[1]))
colnames(trends)<-c('lat','lon','trends','odds','auto')

library(nlme) # this is a package for fitting regressions
# if you do not have it, install it as before

# we are going to cycle over unique lat-lon combinations
# and compute trends now 

for(a in 1:(length(lat_steps)-1)){ # for each latitude step 
	for(b in 1:(length(lon_steps)-1)){ # for each longitude step 
		temp<-dataset[which(dataset$lat==lat_steps[a] & dataset$lon==lon_steps[b]),]$sst	# fetch SST data 
			if( length(which(is.na(temp)==F))>3 & sd(temp,na.rm=T)>0){ # Check that the data is good
			chosen_row<-which(trends$lat==lat_steps[a] & trends$lon==lon_steps[b])
			df1<-as.data.frame(cbind(temp,years))
			df1<-df1[complete.cases(df1),]
			mod1<-gls( temp ~ years ,data=df1) # Fit a standard GLS model 
			tryCatch( mod2<-gls( temp ~ years  ,correlation=corAR1(form=~years),data= df1) ,error=function(e) print('test') )
			# fit a GLS model with temporal autocorrelation; if the fit doesn't work 'test' will be printed to the terminal
			pval<-1
			if(exists('mod2')){
					pval<-anova(mod1,mod2)$p[2] # is the standard or the autocorrelated model preferred?
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
	print(a) # Print the latitude step we are on (1 to 15 steps)
}




library(pracma) # This is a package we need for the regridding of the results
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
# 'Chlor_3' is a reprojected summary of the temporal trends from the preferred model fits

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

library(ggplot2) # package for producing aesthetic plots
# if you do not have this, install it as before 

SST_evolution<-
ggplot(data=Chlor_3, aes(x=x,y=y))+
coord_equal()+
geom_raster(aes(x=x,y=y,fill=z),interpolate=F)+
scale_fill_gradient2(high='red',low='blue',mid='white',limits=limits,na.value='white')+
geom_raster(data = dplyr::filter(bath4, !is.na(z)), 
            aes(x = x, y = y), fill = "gray40")+
geom_point(data=as.data.frame(coordinates(stip_d.new)),
aes(x=lon,y=lat),alpha=.5,size=1.2,shape=point_shape)+
geom_text(data=as.data.frame(coordinates(d.new)),
aes(x=lon,y=lat), label=c(expression(paste('0'^o,'E')),expression(paste('30'^o,'E')),expression(paste('60'^o,'E')),expression(paste('75'^o,'N')),expression(paste('80'^o,'N'))),
angle=c(-30,0,30,30,30))+
labs(fill=expression(paste(Delta, italic('SST'),' y'^-1)) )+
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
# type 'SST_evolution' into the terminal to view the plot
# significant trends at a 0.05 alpha level in SST over time are marked with a circle (standard fit) or triangle (autocorrelated)

# We are going to now perform the same exercise for the remaining variables
# I will desist comments for now 

trends<-cbind(expand.grid(lat_steps,lon_steps),NA,NA)
trends<-as.data.frame(trends)
trends<-cbind(trends,rep(0,dim(trends)[1]))
colnames(trends)<-c('lat','lon','trends','odds','auto')


for(a in 1:(length(lat_steps)-1)){
	for(b in 1:(length(lon_steps)-1)){
		biomass<-dataset[which(dataset$lat==lat_steps[a] & dataset$lon==lon_steps[b]),]$chla
			if( length(which(is.na(biomass)==F))>3 & sd(biomass,na.rm=T)>0){
			chosen_row<-which(trends$lat==lat_steps[a] & trends$lon==lon_steps[b])
			df1<-as.data.frame(cbind(biomass,years))
			df1<-df1[complete.cases(df1),]
			mod1<-gls( biomass ~ years ,data=df1)
			tryCatch( mod2<-gls( biomass ~ years  ,correlation=corAR1(form=~years),data= df1) ,error=function(e) print('test') )
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

Biomass_evolution<-
ggplot(data=Chlor_3, aes(x=x,y=y))+
coord_equal()+
geom_raster(aes(x=x,y=y,fill=z),interpolate=F)+
scale_fill_gradient2(high='red',low='blue',mid='white',limits=limits,na.value='white')+
geom_raster(data = dplyr::filter(bath4, !is.na(z)), 
            aes(x = x, y = y), fill = "gray40")+
geom_point(data=as.data.frame(coordinates(stip_d.new)),
aes(x=lon,y=lat),alpha=.5,size=1.2,shape=point_shape)+
geom_text(data=as.data.frame(coordinates(d.new)),
aes(x=lon,y=lat), label=c(expression(paste('0'^o,'E')),expression(paste('30'^o,'E')),expression(paste('60'^o,'E')),expression(paste('75'^o,'N')),expression(paste('80'^o,'N'))),
angle=c(-30,0,30,30,30))+
labs(fill=expression(paste(Delta, italic('B'),' y'^-1)) )+
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
colnames(trends)<-c('lat','lon','trends','odds','auto')


for(a in 1:(length(lat_steps)-1)){
	for(b in 1:(length(lon_steps)-1)){
		wind<-dataset[which(dataset$lat==lat_steps[a] & dataset$lon==lon_steps[b]),]$storminess
			if( length(which(is.na(wind)==F))>3 & sd(wind,na.rm=T)>0){
			chosen_row<-which(trends$lat==lat_steps[a] & trends$lon==lon_steps[b])
			df1<-as.data.frame(cbind(wind,years))
			df1<-df1[complete.cases(df1),]
			mod1<-gls( wind ~ years ,data=df1)
			tryCatch( mod2<-gls( wind ~ years  ,correlation=corAR1(form=~years),data= df1) ,error=function(e) print('test') )
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

Wind_evolution<-
ggplot(data=Chlor_3, aes(x=x,y=y))+
coord_equal()+
geom_raster(aes(x=x,y=y,fill=z),interpolate=F)+
scale_fill_gradient2(high='red',low='blue',mid='white',limits=limits,na.value='white')+
geom_raster(data = dplyr::filter(bath4, !is.na(z)), 
            aes(x = x, y = y), fill = "gray40")+
geom_point(data=as.data.frame(coordinates(stip_d.new)),
aes(x=lon,y=lat),alpha=.5,size=1.2,shape=point_shape)+
geom_text(data=as.data.frame(coordinates(d.new)),
aes(x=lon,y=lat), label=c(expression(paste('0'^o,'E')),expression(paste('30'^o,'E')),expression(paste('60'^o,'E')),expression(paste('75'^o,'N')),expression(paste('80'^o,'N'))),
angle=c(-30,0,30,30,30))+
labs(fill=expression(paste(Delta, italic('W'),' y'^-1)) )+
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
colnames(trends)<-c('lat','lon','trends','odds','auto')


for(a in 1:(length(lat_steps)-1)){
	for(b in 1:(length(lon_steps)-1)){
		currents<-dataset[which(dataset$lat==lat_steps[a] & dataset$lon==lon_steps[b]),]$velocity
			if( length(which(is.na(currents)==F))>3 & sd(currents,na.rm=T)>0){
			chosen_row<-which(trends$lat==lat_steps[a] & trends$lon==lon_steps[b])
			df1<-as.data.frame(cbind(currents,years))
			df1<-df1[complete.cases(df1),]
			mod1<-gls( currents ~ years ,data=df1)
			tryCatch( mod2<-gls( currents ~ years  ,correlation=corAR1(form=~years),data= df1) ,error=function(e) print('test') )
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

Current_evolution<-
ggplot(data=Chlor_3, aes(x=x,y=y))+
coord_equal()+
geom_raster(aes(x=x,y=y,fill=z),interpolate=F)+
scale_fill_gradient2(high='red',low='blue',mid='white',limits=limits,na.value='white')+
geom_raster(data = dplyr::filter(bath4, !is.na(z)), 
            aes(x = x, y = y), fill = "gray40")+
geom_point(data=as.data.frame(coordinates(stip_d.new)),
aes(x=lon,y=lat),alpha=.5,size=1.2,shape=point_shape)+
geom_text(data=as.data.frame(coordinates(d.new)),
aes(x=lon,y=lat), label=c(expression(paste('0'^o,'E')),expression(paste('30'^o,'E')),expression(paste('60'^o,'E')),expression(paste('75'^o,'N')),expression(paste('80'^o,'N'))),
angle=c(-30,0,30,30,30))+
labs(fill=expression(paste(Delta, italic('C'),' y'^-1)) )+
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

# Now it is time to assemble the final figure 
library(ggpubr)# This is a package for combining plots
# if you do not have it, install it as before


dev.new(width=10,height=10,units='inches') # Make a new plot 
ggarrange(Biomass_evolution,SST_evolution,Wind_evolution,Current_evolution,ncol=2,nrow=2,labels=c('a','b','c','d'),label.x=0.175,label.y=0.8,font.label=list(size=24),
align='hv')
#ggsave(filename='Linear_trends_13_07_2021.pdf',dpi=300)
# If you wish to save the plot, uncomment the above line 



