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

emi_nh3 <- raster('1km_NH3_em_ms10.tif')
emi_nox_1 <- raster('1km_NOx_em_ms2.tif')
emi_nox_2 <- raster('1km_NOx_em_ms34.tif')
emi_nox_3 <- raster('1km_NOx_em_ms7.tif')
emi_nox_4 <-raster('1km_NOx_em_ms8.tif')
emi_pm25_1 <- raster('1km_PM25_em_ms2.tif')
emi_pm25_2 <- raster('1km_pm25_em_ms34.tif')
emi_SO2_1 <- raster('1km_SO2_em_ms1.tif')
emi_SO2_2 <- raster('1km_SO2_em_ms2.tif')
emi_SO2_3 <- raster('1km_so2_em_ms34.tif')
emi_voc_1 <- raster('1km_VOC_em_ms2.tif')
emi_voc_2 <- raster('1km_voc_em_ms34.tif')
emi_voc_3 <- raster('1km_VOC_em_ms5.tif')
emi_voc_4 <- raster('1km_VOC_em_ms6.tif')
emi_voc_5 <- raster('1km_VOC_em_ms7.tif')
emi_voc_6 <- raster('1km_VOC_em_ms8.tif')

#reading output data
con_pm25 <- read.csv('2016_pm25_data.csv')
sum(is.na(con_pm25))
summary(con_pm25)
#obtaining list of all cooridnates across Europe
coord_pm25 <- cbind(con_pm25$SamplingPoint_Longitude, con_pm25$SamplingPoint_Latitude)

#exploratory analysis
#exploratory analysis
hist(con_pm25$AQValue)
ggplot(data = con_pm25,aes(coord_pm25[,1],coord_pm25[,2],colour=AQValue))+geom_point(aes(colour=AQValue))+scale_color_gradientn(colours = terrain.colors(5))+labs(x="Longitude",y="Latitude",title = "Air Concentration of pm25")

#extracting the data for required inputs 
results_nh3_pm25 <- extract(emi_nh3, coord_pm25)
results_nox_1_pm25 <- extract(emi_nox_1, coord_pm25)
results_nox_2_pm25 <- extract(emi_nox_2, coord_pm25)
results_nox_3_pm25 <- extract(emi_nox_3, coord_pm25)
results_nox_4_pm25 <- extract(emi_nox_4, coord_pm25)
results_pm25_1_pm25 <- extract(emi_pm25_1, coord_pm25)
results_pm25_2_pm25 <- extract(emi_pm25_2, coord_pm25)
results_so2_1_pm25 <- extract(emi_SO2_1, coord_pm25)
results_so2_2_pm25 <- extract(emi_SO2_2, coord_pm25)
results_so2_3_pm25 <- extract(emi_SO2_3, coord_pm25)
results_voc_1_pm25 <- extract(emi_voc_1, coord_pm25)
results_voc_2_pm25 <- extract(emi_voc_2, coord_pm25)
results_voc_4_pm25 <- extract(emi_voc_4, coord_pm25)
results_voc_3_pm25 <- extract(emi_voc_3, coord_pm25)
results_voc_5_pm25 <- extract(emi_voc_5, coord_pm25)
results_voc_6_pm25 <- extract(emi_voc_6, coord_pm25)

final_results_nh3_pm25 <- results_nh3_pm25
final_results_nox_pm25 <- results_nox_1_pm25 + results_nox_2_pm25 + results_nox_3_pm25 + results_nox_4_pm25
final_results_pm25_pm25 <- results_pm25_1_pm25 + results_pm25_2_pm25
final_results_so2_pm25 <- results_so2_1_pm25 + results_so2_2_pm25 + results_so2_3_pm25
final_results_voc_pm25 <- results_voc_1_pm25 + results_voc_2_pm25 + results_voc_3_pm25 + results_voc_4_pm25

#sum(is.na(final_results_nh3_pm25))
#sum(is.na(final_results_nox_pm25))
#sum(is.na(final_results_pm25_pm25))
#sum(is.na(final_results_so2_pm25))
#sum(is.na(final_results_voc_pm25))

#replacing NAs with zero
final_results_nh3_pm25[which(is.na(final_results_nh3_pm25))]<-0.0
final_results_nox_pm25[which(is.na(final_results_nox_pm25))]<-0.0
final_results_pm25_pm25[which(is.na(final_results_pm25_pm25))]<-0.0
final_results_so2_pm25[which(is.na(final_results_so2_pm25))]<-0.0
final_results_voc_pm25[which(is.na(final_results_voc_pm25))]<-0.0

#forming dataframe with location coordinates and inp/out data
new_df_pm25 <- data.frame(con_pm25$SamplingPoint_Longitude,con_pm25$SamplingPoint_Latitude,final_results_nh3_pm25,final_results_nox_pm25,final_results_pm25_pm25,final_results_so2_pm25,final_results_voc_pm25,con_pm25$AQValue)

#converting data frame into spatial points data frame
spdf.con_pm25 <- SpatialPointsDataFrame(coords = coord_pm25,data = new_df_pm25[,3:8],proj4string = crs(emi_nh3))
colnames(spdf.con_pm25@data)<-c("Emissions_nh3","Emissions_nox","Emissions_pm25","Emissions_so2","Emissions_voc","AQValue_pm25")

#distance matrix
dmat_pm25<-gw.dist(dp.locat = coordinates(spdf.con_pm25))
#selecting bandwidth for geographically weighted regression
bw_optimal_bisq_pm25 <- bw.gwr(AQValue_pm25~Emissions_nh3+Emissions_nox+Emissions_pm25+Emissions_so2+Emissions_voc,data = spdf.con_pm25,approach = "CV",kernel = "bisquare",adaptive = TRUE,dMat = dmat_pm25)
#forming gwr.basic model
gwr_mod_bisq_pm25<-gwr.basic(AQValue_pm25~Emissions_nh3+Emissions_nox+Emissions_pm25+Emissions_so2+Emissions_voc,data=spdf.con_pm25,bw=bw_optimal_bisq_pm25,kernel = "bisquare",adaptive = TRUE)
#forming gwr.robust model
gwr_mod_bisq_robust_pm25<-gwr.robust(AQValue_pm25~Emissions_nh3+Emissions_nox+Emissions_pm25+Emissions_so2+Emissions_voc,data=spdf.con_pm25,bw=bw_optimal_bisq_pm25,kernel = "bisquare",adaptive = TRUE)

bw_optimal_gauss_pm25 <- bw.gwr(AQValue_pm25~Emissions_nh3+Emissions_nox+Emissions_pm25+Emissions_so2+Emissions_voc,data = spdf.con_pm25,approach = "CV",kernel = "gaussian",adaptive = TRUE,dMat = dmat_pm25)
gwr_mod_gauss_pm25<-gwr.basic(AQValue_pm25~Emissions_nh3+Emissions_nox+Emissions_pm25+Emissions_so2+Emissions_voc,data=spdf.con_pm25,bw=bw_optimal_gauss_pm25,kernel = "gaussian",adaptive = TRUE)
gwr_mod_gauss_robust_pm25<-gwr.robust(AQValue_pm25~Emissions_nh3+Emissions_nox+Emissions_pm25+Emissions_so2+Emissions_voc,data=spdf.con_pm25,bw=bw_optimal_gauss_pm25,kernel = "gaussian",adaptive = TRUE)

bw_optimal_exp_pm25 <- bw.gwr(AQValue_pm25~Emissions_nh3+Emissions_nox+Emissions_pm25+Emissions_so2+Emissions_voc,data = spdf.con_pm25,approach = "CV",kernel = "exponential",adaptive = TRUE,dMat = dmat_pm25)
gwr_mod_exp_pm25<-gwr.basic(AQValue_pm25~Emissions_nh3+Emissions_nox+Emissions_pm25+Emissions_so2+Emissions_voc,data=spdf.con_pm25,bw=bw_optimal_exp_pm25,kernel = "exponential",adaptive = TRUE)
gwr_mod_exp_robust_pm25<-gwr.robust(AQValue_pm25~Emissions_nh3+Emissions_nox+Emissions_pm25+Emissions_so2+Emissions_voc,data=spdf.con_pm25,bw=bw_optimal_exp_pm25,kernel = "exponential",adaptive = TRUE)

gwr_mod_bisq_pm25$GW.diagnostic$gwR2.adj
gwr_mod_bisq_robust_pm25$GW.diagnostic$gwR2.adj
gwr_mod_exp_pm25$GW.diagnostic$gwR2.adj
gwr_mod_exp_robust_pm25$GW.diagnostic$gwR2.adj

spdf.shape<-readOGR("NUTS_RG_10M_2016_4326_LEVL_0.shp")
crs(spdf.shape)
spdf.shape@data$id<-rownames(spdf.shape@data)
spdf.points<-fortify(spdf.shape,region = "id")
str(spdf.points)
countries <- inner_join(spdf.points, spdf.shape@data, by="id")
df_coef_nh3_pm25<-as.data.frame(gwr_mod_bisq_pm25$SDF$Emissions_nh3)
colnames(df_coef_nh3_pm25)<-"nh3"
df_coef_nox_pm25<-as.data.frame(gwr_mod_bisq_pm25$SDF$Emissions_nox)
colnames(df_coef_nox_pm25)<-"NOx"
df_coef_pm25_pm25<-as.data.frame(gwr_mod_bisq_pm25$SDF$Emissions_pm25)
colnames(df_coef_pm25_pm25)<-"pm25"
df_coef_so2_pm25<-as.data.frame(gwr_mod_bisq_pm25$SDF$Emissions_so2)
colnames(df_coef_so2_pm25)<-"so2"
df_coef_voc_pm25<-as.data.frame(gwr_mod_bisq_pm25$SDF$Emissions_voc)
colnames(df_coef_voc_pm25)<-"voc"

summary(df_coef_nh3_pm25)
a_plot_nh3_pm10<-ggplot(data = df_coef_nh3_pm25,aes(coord_pm25[,1],coord_pm25[,2],colour=nh3))+geom_point(aes(colour=nh3))+scale_color_gradient2(low = "green",mid = "white",high = "red")
a_plot_nh3_pm10<-a_plot_nh3_pm10 + geom_polygon(data = countries,colour="black",fill=NA,aes(x=long, y=lat,group=group))
a_plot_nh3_pm10 <- a_plot_nh3_pm10 + labs(x="",y="",title = "Coefficient Estimates of NH3 across Europe")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
a_plot_nh3_pm10

summary(df_coef_nox_pm25)
a_plot_nox_pm10<-ggplot(data = df_coef_nox_pm25,aes(coord_pm25[,1],coord_pm25[,2],colour=NOx))+geom_point(aes(colour=NOx))+scale_color_gradient2(low = "green",mid = "white",high = "red")
a_plot_nox_pm10<-a_plot_nox_pm10 + geom_polygon(data = countries,colour="black",fill=NA,aes(x=long, y=lat,group=group))
a_plot_nox_pm10 <- a_plot_nox_pm10 + labs(x="",y="",title = "Coefficient Estimates of NOx across Europe")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
a_plot_nox_pm10

summary(df_coef_pm25_pm25)
a_plot_pm25<-ggplot(data = df_coef_pm25_pm25,aes(coord_pm25[,1],coord_pm25[,2],colour=pm25))+geom_point(aes(colour=pm25))+scale_color_gradient2(low = "green",mid = "white",high = "red")
a_plot_pm25<-a_plot_pm25 + geom_polygon(data = countries,colour="black",fill=NA,aes(x=long, y=lat,group=group))
a_plot_pm25 <- a_plot_pm25 + labs(x="",y="",title = "Coefficient Estimates of pm25 across Europe")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
a_plot_pm25

summary(df_coef_so2_pm25)
a_plot_so2<-ggplot(data = df_coef_so2_pm25,aes(coord_pm25[,1],coord_pm25[,2],colour=so2))+geom_point(aes(colour=so2))+scale_color_gradient2(low = "green",mid = "white",high = "red")
a_plot_so2<-a_plot_so2 + geom_polygon(data = countries,colour="black",fill=NA,aes(x=long, y=lat,group=group))
a_plot_so2 <- a_plot_so2 + labs(x="",y="",title = "Coefficient Estimates of SO2 across Europe")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
a_plot_so2

summary(df_coef_voc_pm25)
a_plot_voc<-ggplot(data = df_coef_voc_pm25,aes(coord_pm25[,1],coord_pm25[,2],colour=voc))+geom_point(aes(colour=voc))+scale_color_gradient2(low = "green",mid = "white",high = "red")
a_plot_voc<-a_plot_voc + geom_polygon(data = countries,colour="black",fill=NA,aes(x=long, y=lat,group=group))
a_plot_voc <- a_plot_voc + labs(x="",y="",title = "Coefficient Estimates of VOC across Europe")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
a_plot_voc

#-----------------------------------------IRELAND---------------------------------------------

ire_extent<-extent(-10.56,-5.34,51.39,55.43)

emi_nh3_ire <- crop(raster('1km_NH3_em_ms10.tif'),ire_extent)
emi_nox_1_ire <- crop(raster('1km_NOx_em_ms2.tif'),ire_extent)
emi_nox_2_ire <- crop(raster('1km_NOx_em_ms34.tif'),ire_extent)
emi_nox_3_ire <- crop(raster('1km_NOx_em_ms7.tif'),ire_extent)
emi_nox_4_ire <-crop(raster('1km_NOx_em_ms8.tif'),ire_extent)
emi_pm25_1_ire <- crop(raster('1km_pm25_em_ms2.tif'),ire_extent)
emi_pm25_2_ire <- crop(raster('1km_pm25_em_ms34.tif'),ire_extent)
emi_SO2_1_ire <- crop(raster('1km_SO2_em_ms1.tif'),ire_extent)
emi_SO2_2_ire <- crop(raster('1km_SO2_em_ms2.tif'),ire_extent)
emi_SO2_3_ire <- crop(raster('1km_so2_em_ms34.tif'),ire_extent)
emi_voc_1_ire <- crop(raster('1km_VOC_em_ms2.tif'),ire_extent)
emi_voc_2_ire <- crop(raster('1km_voc_em_ms34.tif'),ire_extent)
emi_voc_3_ire <- crop(raster('1km_VOC_em_ms5.tif'),ire_extent)
emi_voc_4_ire <- crop(raster('1km_VOC_em_ms6.tif'),ire_extent)
emi_voc_5_ire <- crop(raster('1km_VOC_em_ms7.tif'),ire_extent)
emi_voc_6_ire <- crop(raster('1km_VOC_em_ms8.tif'),ire_extent)

coord_ire<-coordinates(emi_nh3_ire)

pred_results_nh3_ire <- extract(emi_nh3_ire,coord_ire)
pred_results_nox_1_ire <- extract(emi_nox_1_ire,coord_ire)
pred_results_nox_2_ire <- extract(emi_nox_2_ire,coord_ire)
pred_results_nox_3_ire <- extract(emi_nox_3_ire,coord_ire)
pred_results_nox_4_ire <- extract(emi_nox_4_ire,coord_ire)
pred_results_pm25_1_ire <- extract(emi_pm25_1_ire,coord_ire)
pred_results_pm25_2_ire <- extract(emi_pm25_2_ire,coord_ire)
pred_results_SO2_1_ire <- extract(emi_SO2_1_ire,coord_ire)
pred_results_SO2_2_ire <- extract(emi_SO2_2_ire,coord_ire)
pred_results_SO2_3_ire <- extract(emi_SO2_3_ire,coord_ire)
pred_results_voc_1_ire <- extract(emi_voc_1_ire,coord_ire)
pred_results_voc_2_ire <- extract(emi_voc_2_ire,coord_ire)
pred_results_voc_3_ire <- extract(emi_voc_3_ire,coord_ire)
pred_results_voc_4_ire <- extract(emi_voc_4_ire,coord_ire)
pred_results_voc_5_ire <- extract(emi_voc_5_ire,coord_ire)
pred_results_voc_6_ire <- extract(emi_voc_6_ire,coord_ire)

pred_results_nox_ire<-pred_results_nox_1_ire+pred_results_nox_2_ire+pred_results_nox_3_ire+pred_results_nox_4_ire
pred_results_pm25_ire<-pred_results_pm25_1_ire+pred_results_pm25_2_ire
pred_results_so2_ire<-pred_results_SO2_1_ire+pred_results_SO2_2_ire+pred_results_SO2_3_ire
pred_results_voc_ire<-pred_results_voc_1_ire+pred_results_voc_2_ire+pred_results_voc_3_ire+pred_results_voc_4_ire+pred_results_voc_5_ire+pred_results_voc_6_ire

#reading output data
con <- read.csv('2016_pm25_data.csv')
con_ire_pm25 <- con[con$CountryOrTerritory=="Ireland",]
str(con_ire_pm25)
sum(is.na(con_ire_pm25))
summary(con_ire_pm25)
max_ire_aq_pm25<-max(con_ire_pm25$AQValue)
#obtaining list of all cooridnates across Europe
coord_ire_pm25 <- cbind(con_ire_pm25$SamplingPoint_Longitude, con_ire_pm25$SamplingPoint_Latitude)

#extracting the data for required inputs 
results_nh3_ire <- extract(emi_nh3_ire,coord_ire_pm25)
results_nox_1_ire <- extract(emi_nox_1_ire,coord_ire_pm25)
results_nox_2_ire <- extract(emi_nox_2_ire,coord_ire_pm25)
results_nox_3_ire <- extract(emi_nox_3_ire,coord_ire_pm25)
results_nox_4_ire <- extract(emi_nox_4_ire,coord_ire_pm25)
results_pm25_1_ire <- extract(emi_pm25_1_ire,coord_ire_pm25)
results_pm25_2_ire <- extract(emi_pm25_2_ire,coord_ire_pm25)
results_so2_1_ire <- extract(emi_SO2_1_ire,coord_ire_pm25)
results_so2_2_ire <- extract(emi_SO2_2_ire,coord_ire_pm25)
results_so2_3_ire <- extract(emi_SO2_3_ire,coord_ire_pm25)
results_voc_1_ire <- extract(emi_voc_1_ire,coord_ire_pm25)
results_voc_2_ire <- extract(emi_voc_2_ire,coord_ire_pm25)
results_voc_3_ire <- extract(emi_voc_3_ire,coord_ire_pm25)
results_voc_4_ire <- extract(emi_voc_4_ire,coord_ire_pm25)
results_voc_5_ire <- extract(emi_voc_5_ire,coord_ire_pm25)
results_voc_6_ire <- extract(emi_voc_6_ire,coord_ire_pm25)

#adding inputs from all results and cumulating from all macrosectors
final_results_nh3_ire_pm25 = results_nh3_ire
final_results_nox_ire_pm25 = results_nox_1_ire + results_nox_2_ire + results_nox_3_ire + results_nox_4_ire
final_results_pm25_ire_pm25 = results_pm25_1_ire + results_pm25_2_ire
final_results_so2_ire_pm25 = results_so2_1_ire + results_so2_2_ire + results_so2_3_ire
final_results_voc_ire_pm25 = results_voc_1_ire + results_voc_2_ire + results_voc_3_ire + results_voc_4_ire + results_voc_5_ire + results_voc_6_ire

#replacing NAs with zeroes in final results
final_results_nh3_ire_pm25[which(is.na(final_results_nh3_ire_pm25))]<-0.0
final_results_nox_ire_pm25[which(is.na(final_results_nox_ire_pm25))]<-0.0
final_results_pm25_ire_pm25[which(is.na(final_results_pm25_ire_pm25))]<-0.0
final_results_so2_ire_pm25[which(is.na(final_results_so2_ire_pm25))]<-0.0
final_results_voc_ire_pm25[which(is.na(final_results_voc_ire_pm25))]<-0.0

#forming a new dataframe with coordinates,inputs and outputs
new_df_ire_pm25 <- data.frame(con_ire_pm25$SamplingPoint_Longitude,con_ire_pm25$SamplingPoint_Latitude,final_results_nh3_ire_pm25,final_results_nox_ire_pm25,final_results_pm25_ire_pm25,final_results_so2_ire_pm25,final_results_voc_ire_pm25,con_ire_pm25$AQValue)

#converting data frame into spatial points data frame
spdf.con_ire.pm25 <- SpatialPointsDataFrame(coords = coord_ire_pm25,data = new_df_ire_pm25[,3:8],proj4string = crs(emi_nh3))
colnames(spdf.con_ire.pm25@data)<-c("Emissions_nh3","Emissions_nox","Emissions_pm25","Emissions_so2","Emissions_voc","AQValue_pm25")
spdf.pred_ire.pm25 <- SpatialPointsDataFrame(coords = coord_ire,data = as.data.frame(cbind(pred_results_nh3_ire,pred_results_nox_ire,pred_results_pm25_ire,pred_results_so2_ire,pred_results_voc_ire)),proj4string = crs(emi_nh3))
colnames(spdf.pred_ire.pm25@data)<-c("Emissions_nh3","Emissions_nox","Emissions_pm25","Emissions_so2","Emissions_voc")

dmat_ire_pm25<-gw.dist(dp.locat = coordinates(spdf.con_ire.pm25))
#selecting bandwidth for geographically weighted regression
bw_optimal_gauss_ire_pm25 <- bw.gwr(AQValue_pm25~Emissions_nh3+Emissions_nox+Emissions_pm25+Emissions_so2+Emissions_voc,data = spdf.con_ire.pm25,approach = "CV",kernel = "bisquare",adaptive = TRUE,dMat = dmat_ire_pm25)
#prediction
gwr.pred_ire.pm25<-gwr.predict(AQValue_pm25~Emissions_nh3+Emissions_nox+Emissions_pm25+Emissions_so2+Emissions_voc,data=spdf.con_ire.pm25,predictdata = spdf.pred_ire.pm25,bw=bw_optimal_gauss_ire_pm25,kernel = "bisquare",adaptive = TRUE, dMat1 = gw.dist(dp.locat = coordinates(spdf.con_ire.pm25),rp.locat = coordinates(spdf.pred_ire.pm25)), dMat2 = gw.dist(dp.locat = coordinates(spdf.con_ire.pm25)))
pred_df_ire_pm25<-as.data.frame(gwr.pred_ire.pm25$SDF$prediction)
colnames(pred_df_ire_pm25)<-"Pred_AQValue"
pred_df_ire_pm25$Pred_AQValue[which(pred_df_ire_pm25$Pred_AQValue<0)]=0
pred_df_ire_pm25$Pred_AQValue[which(pred_df_ire_pm25$Pred_AQValue>max_ire_aq_pm25)]=max_ire_aq_pm25
summary(pred_df_ire_pm25)
hist(pred_df_ire_pm25$Pred_AQValue)

#forming gwr.basic model
spdf.shape_ire<-readOGR("Census2011_Admin_Counties_generalised20m.shp")
crs(spdf.shape_ire)
spdf.shape_ire@data$id<-rownames(spdf.shape_ire@data)
spdf.points_ire<-fortify(spdf.shape_ire,region = "id")
counties <- inner_join(spdf.points_ire, spdf.shape_ire@data, by="id")
tm65 <- "+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=1.000035 +x_0=200000 +y_0=250000 +a=6377340.189 +b=6356034.447938534 +units=m +no_defs "
newlonglat_ire<-project(cbind(counties$long,counties$lat),proj = tm65,inverse = TRUE)
counties$long<-newlonglat_ire[,1]
counties$lat<-newlonglat_ire[,2]
summary(pred_df_ire_pm25)
k_plot<-ggplot(data=pred_df_ire_pm25,aes(x=coord_ire[,1],y=coord_ire[,2]))
k_plot<-k_plot+geom_point(aes(colour=Pred_AQValue))+scale_color_gradient2(midpoint = 3,low = "green",mid = "white",high = "red")+geom_polygon(data = counties,colour="black",fill=NA,aes(x=long, y=lat, group=group))
k_plot<-k_plot + labs(x="",y="",title = "Predicted Conc. of pm25 in Ireland")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
k_plot

#-----------------------------------------ITALY---------------------------------------------

it_extent<-extent(6.7499552751,18.4802470232, 36.619987291, 47.1153931748)

emi_nh3_it <- crop(raster('1km_NH3_em_ms10.tif'),it_extent)
emi_nox_1_it <- crop(raster('1km_NOx_em_ms2.tif'),it_extent)
emi_nox_2_it <- crop(raster('1km_NOx_em_ms34.tif'),it_extent)
emi_nox_3_it <- crop(raster('1km_NOx_em_ms7.tif'),it_extent)
emi_nox_4_it <-crop(raster('1km_NOx_em_ms8.tif'),it_extent)
emi_pm25_1_it <- crop(raster('1km_PM25_em_ms2.tif'),it_extent)
emi_pm25_2_it <- crop(raster('1km_pm25_em_ms34.tif'),it_extent)
emi_SO2_1_it <- crop(raster('1km_SO2_em_ms1.tif'),it_extent)
emi_SO2_2_it <- crop(raster('1km_SO2_em_ms2.tif'),it_extent)
emi_SO2_3_it <- crop(raster('1km_so2_em_ms34.tif'),it_extent)
emi_voc_1_it <- crop(raster('1km_VOC_em_ms2.tif'),it_extent)
emi_voc_2_it <- crop(raster('1km_voc_em_ms34.tif'),it_extent)
emi_voc_3_it <- crop(raster('1km_VOC_em_ms5.tif'),it_extent)
emi_voc_4_it <- crop(raster('1km_VOC_em_ms6.tif'),it_extent)
emi_voc_5_it <- crop(raster('1km_VOC_em_ms7.tif'),it_extent)
emi_voc_6_it <- crop(raster('1km_VOC_em_ms8.tif'),it_extent)

coord_it<-coordinates(emi_nh3_it)

pred_results_nh3_it <- extract(emi_nh3_it,coord_it)
pred_results_nox_1_it <- extract(emi_nox_1_it,coord_it)
pred_results_nox_2_it <- extract(emi_nox_2_it,coord_it)
pred_results_nox_3_it <- extract(emi_nox_3_it,coord_it)
pred_results_nox_4_it <- extract(emi_nox_4_it,coord_it)
pred_results_pm25_1_it <- extract(emi_pm25_1_it,coord_it)
pred_results_pm25_2_it <- extract(emi_pm25_2_it,coord_it)
pred_results_SO2_1_it <- extract(emi_SO2_1_it,coord_it)
pred_results_SO2_2_it <- extract(emi_SO2_2_it,coord_it)
pred_results_SO2_3_it <- extract(emi_SO2_3_it,coord_it)
pred_results_voc_1_it <- extract(emi_voc_1_it,coord_it)
pred_results_voc_2_it <- extract(emi_voc_2_it,coord_it)
pred_results_voc_3_it <- extract(emi_voc_3_it,coord_it)
pred_results_voc_4_it <- extract(emi_voc_4_it,coord_it)
pred_results_voc_5_it <- extract(emi_voc_5_it,coord_it)
pred_results_voc_6_it <- extract(emi_voc_6_it,coord_it)

pred_results_nox_it<-pred_results_nox_1_it+pred_results_nox_2_it+pred_results_nox_3_it+pred_results_nox_4_it
pred_results_pm25_it<-pred_results_pm25_1_it+pred_results_pm25_2_it
pred_results_so2_it<-pred_results_SO2_1_it+pred_results_SO2_2_it+pred_results_SO2_3_it
pred_results_voc_it<-pred_results_voc_1_it+pred_results_voc_2_it+pred_results_voc_3_it+pred_results_voc_4_it+pred_results_voc_5_it+pred_results_voc_6_it

#reading output data
con <- read.csv('2016_pm25_data.csv')
con_it_pm25 <- con[con$CountryOrTerritory=="Italy",]
str(con_it_pm25)
sum(is.na(con_it_pm25))
summary(con_it_pm25)
max_it_aq_pm25<-max(con_it_pm25$AQValue)
#obtaining list of all cooridnates across Europe
coord_it_pm25 <- cbind(con_it_pm25$SamplingPoint_Longitude, con_it_pm25$SamplingPoint_Latitude)

#extracting the data for required inputs 
results_nh3_it <- extract(emi_nh3_it,coord_it_pm25)
results_nox_1_it <- extract(emi_nox_1_it,coord_it_pm25)
results_nox_2_it <- extract(emi_nox_2_it,coord_it_pm25)
results_nox_3_it <- extract(emi_nox_3_it,coord_it_pm25)
results_nox_4_it <- extract(emi_nox_4_it,coord_it_pm25)
results_pm25_1_it <- extract(emi_pm25_1_it,coord_it_pm25)
results_pm25_2_it <- extract(emi_pm25_2_it,coord_it_pm25)
results_so2_1_it <- extract(emi_SO2_1_it,coord_it_pm25)
results_so2_2_it <- extract(emi_SO2_2_it,coord_it_pm25)
results_so2_3_it <- extract(emi_SO2_3_it,coord_it_pm25)
results_voc_1_it <- extract(emi_voc_1_it,coord_it_pm25)
results_voc_2_it <- extract(emi_voc_2_it,coord_it_pm25)
results_voc_3_it <- extract(emi_voc_3_it,coord_it_pm25)
results_voc_4_it <- extract(emi_voc_4_it,coord_it_pm25)
results_voc_5_it <- extract(emi_voc_5_it,coord_it_pm25)
results_voc_6_it <- extract(emi_voc_6_it,coord_it_pm25)

#adding inputs from all results and cumulating from all macrosectors
final_results_nh3_it_pm25 = results_nh3_it
final_results_nox_it_pm25 = results_nox_1_it + results_nox_2_it + results_nox_3_it + results_nox_4_it
final_results_pm25_it_pm25 = results_pm25_1_it + results_pm25_2_it
final_results_so2_it_pm25 = results_so2_1_it + results_so2_2_it + results_so2_3_it
final_results_voc_it_pm25 = results_voc_1_it + results_voc_2_it + results_voc_3_it + results_voc_4_it + results_voc_5_it + results_voc_6_it

#replacing NAs with zeroes in final results
final_results_nh3_it_pm25[which(is.na(final_results_nh3_it_pm25))]<-0.0
final_results_nox_it_pm25[which(is.na(final_results_nox_it_pm25))]<-0.0
final_results_pm25_it_pm25[which(is.na(final_results_pm25_it_pm25))]<-0.0
final_results_so2_it_pm25[which(is.na(final_results_so2_it_pm25))]<-0.0
final_results_voc_it_pm25[which(is.na(final_results_voc_it_pm25))]<-0.0

#forming a new dataframe with coordinates,inputs and outputs
new_df_it_pm25 <- data.frame(con_it_pm25$SamplingPoint_Longitude,con_it_pm25$SamplingPoint_Latitude,final_results_nh3_it_pm25,final_results_nox_it_pm25,final_results_pm25_it_pm25,final_results_so2_it_pm25,final_results_voc_it_pm25,con_it_pm25$AQValue)

#converting data frame into spatial points data frame
spdf.con_it.pm25 <- SpatialPointsDataFrame(coords = coord_it_pm25,data = new_df_it_pm25[,3:8],proj4string = crs(emi_nh3))
colnames(spdf.con_it.pm25@data)<-c("Emissions_nh3","Emissions_nox","Emissions_pm25","Emissions_so2","Emissions_voc","AQValue_pm25")
spdf.pred_it.pm25 <- SpatialPointsDataFrame(coords = coord_it,data = as.data.frame(cbind(pred_results_nh3_it,pred_results_nox_it,pred_results_pm25_it,pred_results_so2_it,pred_results_voc_it)),proj4string = crs(emi_nh3))
colnames(spdf.pred_it.pm25@data)<-c("Emissions_nh3","Emissions_nox","Emissions_pm25","Emissions_so2","Emissions_voc")

dmat_it_pm25<-gw.dist(dp.locat = coordinates(spdf.con_it.pm25))
#selecting bandwidth for geographically weighted regression
bw_optimal_gauss_it_pm25 <- bw.gwr(AQValue_pm25~Emissions_nh3+Emissions_nox+Emissions_pm25+Emissions_so2+Emissions_voc,data = spdf.con_it.pm25,approach = "CV",kernel = "bisquare",adaptive = TRUE,dMat = dmat_it_pm25)
#prediction
gwr.pred_it.pm25<-gwr.predict(AQValue_pm25~Emissions_nh3+Emissions_nox+Emissions_pm25+Emissions_so2+Emissions_voc,data=spdf.con_it.pm25,predictdata = spdf.pred_it.pm25,bw=bw_optimal_gauss_it_pm25,kernel = "bisquare",adaptive = TRUE, dMat1 = gw.dist(dp.locat = coordinates(spdf.con_it.pm25),rp.locat = coordinates(spdf.pred_it.pm25)), dMat2 = gw.dist(dp.locat = coordinates(spdf.con_it.pm25)))
pred_df_it_pm25<-as.data.frame(gwr.pred_it.pm25$SDF$prediction)
colnames(pred_df_it_pm25)<-"Pred_AQValue"
pred_df_it_pm25$Pred_AQValue[which(pred_df_it_pm25$Pred_AQValue<0)]=0
pred_df_it_pm25$Pred_AQValue[which(pred_df_it_pm25$Pred_AQValue>max_it_aq_pm25)]=max_it_aq_pm25
summary(pred_df_it_pm25)
hist(pred_df_it_pm25$Pred_AQValue)

#forming gwr.basic model
spdf.shape_it<-readOGR("ITA_adm0.shp")
crs(spdf.shape_it)
spdf.shape_it@data$id<-rownames(spdf.shape_it@data)
spdf.points_it<-fortify(spdf.shape_it,region = "id")
states <- inner_join(spdf.points_it, spdf.shape_it@data, by="id")

summary(pred_df_it_pm25)
l_plot<-ggplot(data=pred_df_it_pm25,aes(x=coord_it[,1],y=coord_it[,2]))
l_plot<-l_plot+geom_point(aes(colour=Pred_AQValue))+scale_color_gradient2(midpoint = 10,low = "green",mid = "white",high = "red")+geom_polygon(data = states,colour="black",fill=NA,aes(x=long, y=lat, group=group))
l_plot<-l_plot + labs(x="",y="",title = "Predicted Conc. of pm25 in Italy")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
l_plot

