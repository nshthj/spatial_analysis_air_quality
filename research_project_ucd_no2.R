#install.packages("rgdal")
#install.packages("GWmodel")
#install.packages("rgeos")
library(raster)
library(dplyr)
library(sp)
library(GWmodel)
library(RColorBrewer)
library(ggplot2)
library(rgeos)
library(proj4)
library(maptools)
library(ggmap)
library(rgdal)
#-------code for making a model with the entire dataset-------

#input files reading rasters
emi_1 <- raster('1km_NOx_em_ms2.tif')
emi_2 <- raster('1km_NOx_em_ms34.tif')
emi_3 <- raster('1km_NOx_em_ms7.tif')
emi_4 <- raster('1km_NOx_em_ms8.tif')

#reading output data
con <- read.csv('2016_no2_data.csv')
summary(con)

#obtaining list of all cooridnates across Europe
coord <- cbind(con$SamplingPoint_Longitude, con$SamplingPoint_Latitude)

#obtaining the spatialpoints for the coordinates only present in the output
#spatialPoints_1 <- SpatialPoints(coords = coord, proj4string = crs(emi_1))

#exploratory analysis
hist(con$AQValue)
ggplot(data = con,aes(coord[,1],coord[,2],colour=AQValue))+geom_point(aes(colour=AQValue))+scale_color_gradientn(colours = terrain.colors(5))+labs(x="Longitude",y="Latitude",title = "Air Concentration of NO2")

#extracting data for the required data points
results_1 <- extract(emi_1, coord)
results_2 <- extract(emi_2, coord)
results_3 <- extract(emi_3, coord)
results_4 <- extract(emi_4, coord)

#adding inputs from all results and cumulating from all macrosectors
final_results = results_1+results_2+results_3+results_4
summary(final_results)
sum(is.na(final_results))
#replacing NAs with zeroes in final results
final_results[which(is.na(final_results))]<-0.0

#forming a new dataframe with coordinates,inputs and outputs
new_df <- data.frame(con$SamplingPoint_Longitude,con$SamplingPoint_Latitude,final_results,con$AQValue)
summary(new_df)
#converting data frame into spatial points data frame
spdf.con <- SpatialPointsDataFrame(coords = coord,data = new_df[,3:4],proj4string = crs(emi_2))
colnames(spdf.con@data)<-c("Emissions_nox","AQValue_no2")
#distance matrix
dmat<-gw.dist(dp.locat = coordinates(spdf.con))
#selecting bandwidth for geographically weighted regression
bw_optimal_bisq <- bw.gwr(AQValue_no2~Emissions_nox,data = spdf.con,approach = "CV",kernel = "bisquare",adaptive = TRUE,dMat = dmat)

#running model for bisqaure, gaussian and exponential for comparision and plotting the coefficients of the best one
#forming gwr.basic model
gwr_mod_bisq<-gwr.basic(AQValue_no2~Emissions_nox,data=spdf.con,bw=bw_optimal_bisq,kernel = "bisquare",adaptive = TRUE)
#forming gwr.robust model
gwr_mod_bisq_robust<-gwr.robust(AQValue_no2~Emissions_nox,data=spdf.con,bw=bw_optimal_bisq,kernel = "bisquare",adaptive = TRUE)

#calculating bandwidth for gaussian model
bw_optimal_gauss <- bw.gwr(AQValue_no2~Emissions_nox,data = spdf.con,approach = "CV",kernel = "gaussian",adaptive = TRUE,dMat = dmat)
gwr_mod_gauss<-gwr.basic(AQValue_no2~Emissions_nox,data=spdf.con,bw=bw_optimal_gauss,kernel = "gaussian",adaptive = TRUE)
gwr_mod_gauss_robust<-gwr.robust(AQValue_no2~Emissions_nox,data=spdf.con,bw=bw_optimal_gauss,kernel = "gaussian",adaptive = TRUE)

#calculating bandwidth for exponential model
bw_optimal_exp <- bw.gwr(AQValue_no2~Emissions_nox,data = spdf.con,approach = "CV",kernel = "exponential",adaptive = TRUE,dMat = dmat)
gwr_mod_exp<-gwr.basic(AQValue_no2~Emissions_nox,data=spdf.con,bw=bw_optimal_exp,kernel = "exponential",adaptive = TRUE)
gwr_mod_exp_robust<-gwr.robust(AQValue_no2~Emissions_nox,data=spdf.con,bw=bw_optimal_exp,kernel = "exponential",adaptive = TRUE)

gwr_mod_bisq$GW.diagnostic$gwR2.adj
gwr_mod_bisq_robust$GW.diagnostic$gwR2.adj
gwr_mod_gauss$GW.diagnostic$gwR2.adj
gwr_mod_gauss_robust$GW.diagnostic$gwR2.adj
gwr_mod_exp$GW.diagnostic$gwR2.adj
gwr_mod_exp_robust$GW.diagnostic$gwR2.adj

spdf.shape<-readOGR("NUTS_RG_10M_2016_4326_LEVL_0.shp")
crs(spdf.shape)
spdf.shape@data$id<-rownames(spdf.shape@data)
spdf.points<-fortify(spdf.shape,region = "id")
str(spdf.points)
countries <- inner_join(spdf.points, spdf.shape@data, by="id")
df_coef<-as.data.frame(gwr_mod_bisq$SDF$Emissions_nox)
colnames(df_coef)<-"NOx"
a_plot<-ggplot(data = df_coef,aes(coord[,1],coord[,2],colour=NOx))+geom_point(aes(colour=NOx))+scale_color_gradient2(midpoint = 5, low = "green",mid = "white",high = "red")
a_plot<-a_plot + geom_polygon(data = countries,colour="black",fill=NA,aes(x=long, y=lat,group=group))
a_plot <- a_plot + labs(x="",y="",title = "Coefficient Estimates of NOx across Europe")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
a_plot
#-----------------------------------------IRELAND-----------------------------------------------

ire_extent<-extent(-10.56,-5.34,51.39,55.43)
emi_1_ire<-crop(emi_1,ire_extent)
emi_2_ire<-crop(emi_2,ire_extent)
emi_3_ire<-crop(emi_3,ire_extent)
emi_4_ire<-crop(emi_4,ire_extent)

coord_ire <- coordinates(emi_1_ire)

pred_results_2_ire <- extract(emi_2_ire,coord_ire)
pred_results_1_ire <- extract(emi_1_ire,coord_ire)
pred_results_3_ire <- extract(emi_3_ire,coord_ire)
pred_results_4_ire <- extract(emi_4_ire,coord_ire)

pred_results_ire <- pred_results_1_ire + pred_results_2_ire + pred_results_3_ire + pred_results_4_ire
sum(is.na(pred_results_ire))
summary(pred_results_ire)

#reading output data
con_ire <- con[con$CountryOrTerritory=="Ireland",]
sum(is.na(con_ire))
summary(con_ire)
max_ire_aq<-max(con_ire$AQValue)
#obtaining list of all cooridnates across Europe
coord_ire_1 <- cbind(con_ire$SamplingPoint_Longitude, con_ire$SamplingPoint_Latitude)

#extracting the data for required inputs 
results_1 <- extract(emi_1_ire, coord_ire_1)
results_2 <- extract(emi_2_ire, coord_ire_1)
results_3 <- extract(emi_3_ire, coord_ire_1)
results_4 <- extract(emi_4_ire, coord_ire_1)

#adding inputs from all results and cumulating from all macrosectors
final_results_ire = results_1+results_2+results_3+results_4
summary(final_results_ire)
sum(is.na(final_results))
#replacing NAs with zeroes in final results
final_results[which(is.na(final_results))]<-0.0

#forming a new dataframe with coordinates,inputs and outputs
new_df_ire <- data.frame(con_ire$SamplingPoint_Longitude,con_ire$SamplingPoint_Latitude,final_results_ire,con_ire$AQValue)
summary(new_df_ire)
#converting data frame into spatial points data frame
spdf.con_ire <- SpatialPointsDataFrame(coords = coord_ire_1,data = new_df_ire[,3:4],proj4string = crs(emi_2))
colnames(spdf.con_ire@data)<-c("Emissions_nox","AQValue_no2")
spdf.pred_ire <- SpatialPointsDataFrame(coords = coord_ire,data = as.data.frame(pred_results_ire),proj4string = crs(emi_1))
colnames(spdf.pred_ire@data)<-c("Emissions_nox")

#distance matrix
dmat_ire<-gw.dist(dp.locat = coordinates(spdf.con_ire))
#selecting bandwidth for geographically weighted regression
bw_optimal_bisq_ire <- bw.gwr(AQValue_no2~Emissions_nox,data = spdf.con_ire,approach = "CV",kernel = "bisquare",adaptive = TRUE,dMat = dmat_ire)
#predicted results
gwr.pred_ire<-gwr.predict(AQValue_no2~Emissions_nox,data=spdf.con_ire,predictdata = spdf.pred_ire,bw=bw_optimal_bisq_ire,kernel = "bisquare",adaptive = TRUE, dMat1 = gw.dist(dp.locat = coordinates(spdf.con_ire),rp.locat = coordinates(spdf.pred_ire)), dMat2 = gw.dist(dp.locat = coordinates(spdf.con_ire)))
pred_df_ire<-as.data.frame(gwr.pred_ire$SDF$prediction)
colnames(pred_df_ire)<-"Pred_AQValue"
pred_df_ire$Pred_AQValue[which(pred_df_ire$Pred_AQValue<0)]=0
pred_df_ire$Pred_AQValue[which(pred_df_ire$Pred_AQValue>max_ire_aq)]=max_ire_aq
summary(pred_df_ire)
hist(pred_df_ire$Pred_AQValue)

#plotting the predicted results for Ireland
spdf.shape_ire<-readOGR("Census2011_Admin_Counties_generalised20m.shp")
crs(spdf.shape_ire)
spdf.shape_ire@data$id<-rownames(spdf.shape_ire@data)
spdf.points_ire<-fortify(spdf.shape_ire,region = "id")
counties <- inner_join(spdf.points_ire, spdf.shape_ire@data, by="id")
tm65 <- "+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=1.000035 +x_0=200000 +y_0=250000 +a=6377340.189 +b=6356034.447938534 +units=m +no_defs "
newlonglat_ire<-project(cbind(counties$long,counties$lat),proj = tm65,inverse = TRUE)
counties$long<-newlonglat_ire[,1]
counties$lat<-newlonglat_ire[,2]

b_plot<-ggplot(data=pred_df_ire,aes(x=coord_ire[,1],y=coord_ire[,2]))
b_plot<-b_plot+geom_point(aes(colour=Pred_AQValue))+scale_color_gradient2(midpoint = 10,low = "green",mid = "white",high = "red")+geom_polygon(data = counties,colour="black",fill=NA,aes(x=long, y=lat, group=group))
b_plot<-b_plot + labs(x="",y="",title = "Predicted AQValue (NO2) in Ireland")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
b_plot

#-------------------------------------------ITALY---------------------------------------------

it_extent<-extent(6.7499552751,18.4802470232, 36.619987291, 47.1153931748)
emi_1_it<-crop(emi_1,it_extent)
emi_2_it<-crop(emi_2,it_extent)
emi_3_it<-crop(emi_3,it_extent)
emi_4_it<-crop(emi_4,it_extent)

coord_it <- coordinates(emi_1_it)

pred_results_1_it <- extract(emi_1_it,coord_it)
pred_results_2_it <- extract(emi_2_it,coord_it)
pred_results_3_it <- extract(emi_3_it,coord_it)
pred_results_4_it <- extract(emi_4_it,coord_it)

pred_results_it <- pred_results_1_it + pred_results_2_it + pred_results_3_it + pred_results_4_it
sum(is.na(pred_results_it))
summary(pred_results_it)

#reading output data
con_it <- con[con$CountryOrTerritory=="Italy",]
sum(is.na(con_it))
summary(con_it)
max_it_aq<-max(con_it$AQValue)
#obtaining list of all cooridnates across Europe
coord_it_1 <- cbind(con_it$SamplingPoint_Longitude, con_it$SamplingPoint_Latitude)

#extracting the data for required inputs 
results_1_it <- extract(emi_1_it, coord_it_1)
results_2_it <- extract(emi_2_it, coord_it_1)
results_3_it <- extract(emi_3_it, coord_it_1)
results_4_it <- extract(emi_4_it, coord_it_1)

#adding inputs from all results and cumulating from all macrosectors
final_results_it = results_1_it+results_2_it+results_3_it+results_4_it
summary(final_results_it)
sum(is.na(final_results))
#replacing NAs with zeroes in final results
final_results[which(is.na(final_results))]<-0.0

#forming a new dataframe with coordinates,inputs and outputs
new_df_it <- data.frame(con_it$SamplingPoint_Longitude,con_it$SamplingPoint_Latitude,final_results_it,con_it$AQValue)
summary(new_df_it)
#converting data frame into spatial points data frame
spdf.con_it <- SpatialPointsDataFrame(coords = coord_it_1,data = new_df_it[,3:4],proj4string = crs(emi_2))
colnames(spdf.con_it@data)<-c("Emissions_nox","AQValue_no2")
spdf.pred_it <- SpatialPointsDataFrame(coords = coord_it,data = as.data.frame(pred_results_it),proj4string = crs(emi_1))
colnames(spdf.pred_it@data)<-c("Emissions_nox")

#distance matrix
dmat_it<-gw.dist(dp.locat = coordinates(spdf.con_it))
#selecting bandwidth for geographically weighted regression
bw_optimal_bisq_it <- bw.gwr(AQValue_no2~Emissions_nox,data = spdf.con_it,approach = "CV",kernel = "bisquare",adaptive = TRUE,dMat = dmat_it)
#predicted results
gwr.pred_it<-gwr.predict(AQValue_no2~Emissions_nox,data=spdf.con_it,predictdata = spdf.pred_it,bw=bw_optimal_bisq_it,kernel = "bisquare",adaptive = TRUE, dMat1 = gw.dist(dp.locat = coordinates(spdf.con_it),rp.locat = coordinates(spdf.pred_it)), dMat2 = gw.dist(dp.locat = coordinates(spdf.con_it)))
pred_df_it<-as.data.frame(gwr.pred_it$SDF$prediction)
colnames(pred_df_it)<-"Pred_AQValue"
pred_df_it$Pred_AQValue[which(pred_df_it$Pred_AQValue<0)]=0
pred_df_it$Pred_AQValue[which(pred_df_it$Pred_AQValue>max_it_aq)]=max_it_aq
summary(pred_df_it)
hist(pred_df_it$Pred_AQValue)

#plotting the predicted results for Ireland
spdf.shape_it<-readOGR("ITA_adm0.shp")
crs(spdf.shape_it)
spdf.shape_it@data$id<-rownames(spdf.shape_it@data)
spdf.points_it<-fortify(spdf.shape_it,region = "id")
states <- inner_join(spdf.points_it, spdf.shape_it@data, by="id")
#tm65 <- "+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=1.000035 +x_0=200000 +y_0=250000 +a=6377340.189 +b=6356034.447938534 +units=m +no_defs "
#newlonglat_it<-project(cbind(counties$long,counties$lat),proj = tm65,inverse = TRUE)
#states$long<-newlonglat_it[,1]
#states$lat<-newlonglat_it[,2]

c_plot<-ggplot(data=pred_df_it,aes(x=coord_it[,1],y=coord_it[,2]))
c_plot<-c_plot+geom_point(aes(colour=Pred_AQValue))+scale_color_gradient2(midpoint = 15,low = "green",mid = "white",high = "red")+geom_polygon(data = states,colour="black",fill=NA,aes(x=long, y=lat, group=group))
c_plot<-c_plot + labs(x="",y="",title = "Predicted AQValue (NO2) in Italy")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
c_plot
