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

#input files reading rasters
emi_nox_1 <- raster('1km_NOx_em_ms2.tif')
emi_nox_2 <- raster('1km_NOx_em_ms34.tif')
emi_nox_3 <- raster('1km_NOx_em_ms7.tif')
emi_nox_4 <- raster('1km_NOx_em_ms8.tif')
emi_voc_1 <- raster('1km_VOC_em_ms2.tif')
emi_voc_2 <- raster('1km_voc_em_ms34.tif')
emi_voc_3 <- raster('1km_VOC_em_ms5.tif')
emi_voc_4 <- raster('1km_VOC_em_ms6.tif')
emi_voc_5 <- raster('1km_VOC_em_ms7.tif')
emi_voc_6 <- raster('1km_VOC_em_ms8.tif')

#reading output data
con_somo <- read.csv('2016_somo35_data.csv')
str(con_somo)
sum(is.na(con_somo))
summary(con_somo)


#obtaining list of all cooridnates across Europe
coord_somo <- cbind(con_somo$SamplingPoint_Longitude, con_somo$SamplingPoint_Latitude)

#exploratory analysis
hist(con$AQValue)
ggplot(data = con_somo,aes(coord_somo[,1],coord_somo[,2],colour=AQValue))+geom_point(aes(colour=AQValue))+scale_color_gradientn(colours = terrain.colors(5))+labs(x="Longitude",y="Latitude",title = "Air Conc. of SOMO35")

#extracting the data for required inputs 
results_nox_1 <- extract(emi_nox_1, coord_somo)
results_nox_2 <- extract(emi_nox_2, coord_somo)
results_nox_3 <- extract(emi_nox_3, coord_somo)
results_nox_4 <- extract(emi_nox_4, coord_somo)
results_voc_1 <- extract(emi_voc_1, coord_somo)
results_voc_2 <- extract(emi_voc_2, coord_somo)
results_voc_3 <- extract(emi_voc_3, coord_somo)
results_voc_4 <- extract(emi_voc_4, coord_somo)
results_voc_5 <- extract(emi_voc_5, coord_somo)
results_voc_6 <- extract(emi_voc_6, coord_somo)

#adding inputs from all results and cumulating from all macrosectors
final_results_nox = results_nox_1 + results_nox_2 + results_nox_3 + results_nox_4
final_results_voc = results_voc_1 + results_voc_2 + results_voc_3 + results_voc_4 + results_voc_5 + results_voc_6

#replacing NAs with zeroes in final results
final_results_nox[which(is.na(final_results_nox))]<-0.0
final_results_voc[which(is.na(final_results_voc))]<-0.0

#forming a new dataframe with coordinates,inputs and outputs
new_df_somo <- data.frame(con_somo$SamplingPoint_Longitude,con_somo$SamplingPoint_Latitude,final_results_nox,final_results_voc,con_somo$AQValue)

#converting data frame into spatial points data frame
spdf.con_somo <- SpatialPointsDataFrame(coords = coord_somo,data = new_df_somo[,3:5],proj4string = crs(emi_nox_1))
colnames(spdf.con_somo@data)<-c("Emissions_nox","Emissions_voc","AQValue_somo35")

dmat_somo<-gw.dist(dp.locat = coordinates(spdf.con_somo))
#selecting bandwidth for geographically weighted regression
bw_optimal_bisq_somo <- bw.gwr(AQValue_somo35~Emissions_nox+Emissions_voc,data = spdf.con_somo,approach = "CV",kernel = "bisquare",adaptive = TRUE,dMat = dmat_somo)
#forming gwr.basic model
gwr_mod_bisq_somo<-gwr.basic(AQValue_somo35~Emissions_nox+Emissions_voc,data=spdf.con_somo,bw=bw_optimal_bisq_somo,kernel = "bisquare",adaptive = TRUE)
#forming gwr.robust model
gwr_mod_bisq_robust_somo<-gwr.robust(AQValue_somo35~Emissions_nox+Emissions_voc,data=spdf.con_somo,bw=bw_optimal_bisq_somo,kernel = "bisquare",adaptive = TRUE)

bw_optimal_gauss_somo <- bw.gwr(AQValue_somo35~Emissions_nox+Emissions_voc,data = spdf.con_somo,approach = "CV",kernel = "gaussian",adaptive = TRUE,dMat = dmat_somo)
gwr_mod_gauss_somo<-gwr.basic(AQValue_somo35~Emissions_nox+Emissions_voc,data=spdf.con_somo,bw=bw_optimal_gauss_somo,kernel = "gaussian",adaptive = TRUE)
gwr_mod_gauss_robust_somo<-gwr.robust(AQValue_somo35~Emissions_nox+Emissions_voc,data=spdf.con_somo,bw=bw_optimal_gauss_somo,kernel = "gaussian",adaptive = TRUE)

bw_optimal_exp_somo <- bw.gwr(AQValue_somo35~Emissions_nox+Emissions_voc,data = spdf.con_somo,approach = "CV",kernel = "exponential",adaptive = TRUE,dMat = dmat_somo)
gwr_mod_exp_somo<-gwr.basic(AQValue_somo35~Emissions_nox+Emissions_voc,data=spdf.con_somo,bw=bw_optimal_exp_somo,kernel = "exponential",adaptive = TRUE)
gwr_mod_exp_robust_somo<-gwr.robust(AQValue_somo35~Emissions_nox+Emissions_voc,data=spdf.con_somo,bw=bw_optimal_exp_somo,kernel = "exponential",adaptive = TRUE)

gwr_mod_bisq_somo$GW.diagnostic$gwR2.adj
gwr_mod_exp_somo$GW.diagnostic$gwR2.adj
gwr_mod_exp_robust_somo$GW.diagnostic$gwR2.adj

spdf.shape_somo<-readOGR("NUTS_RG_10M_2016_4326_LEVL_0.shp")
crs(spdf.shape)
spdf.shape_somo@data$id<-rownames(spdf.shape_somo@data)
spdf.points_somo<-fortify(spdf.shape_somo,region = "id")
str(spdf.points_somo)
countries <- inner_join(spdf.points, spdf.shape@data, by="id")
df_coef_nox<-as.data.frame(gwr_mod_bisq_somo$SDF$Emissions_nox)
df_coef_voc<-as.data.frame(gwr_mod_bisq_somo$SDF$Emissions_voc)
colnames(df_coef_nox)<-"NOx"
colnames(df_coef_voc)<-"VOC"
d_plot<-ggplot(data = df_coef_nox,aes(coord_somo[,1],coord_somo[,2],colour=NOx))+geom_point(aes(colour=NOx))+scale_color_gradient2(low = "green",mid = "white",high = "red")
d_plot<-d_plot + geom_polygon(data = countries,colour="black",fill=NA,aes(x=long, y=lat,group=group))
d_plot <- d_plot + labs(x="",y="",title = "Coefficient Estimates of NOx across Europe")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
d_plot
length(df_coef_voc$VOC)
length(which(df_coef_voc$VOC< -100))
gwr_mod_bisq_somo
e_plot<-ggplot(data = df_coef_voc,aes(coord_somo[,1],coord_somo[,2],colour=VOC))+geom_point(aes(colour=VOC))+scale_color_gradient2(midpoint=-1000,low="green",mid = "white",high = "red")
e_plot<-e_plot + geom_polygon(data = countries,colour="black",fill=NA,aes(x=long, y=lat,group=group))
e_plot <- e_plot + labs(x="",y="",title = "Coefficient Estimates of VOC across Europe")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
e_plot

#----------------------------------------------IRELAND-------------------------------------------

ire_extent<-extent(-10.56,-5.34,51.39,55.43)
emi_nox_1_ire<-crop(emi_nox_1,ire_extent)
emi_nox_2_ire<-crop(emi_nox_2,ire_extent)
emi_nox_3_ire<-crop(emi_nox_3,ire_extent)
emi_nox_4_ire<-crop(emi_nox_4,ire_extent)
emi_voc_1_ire<-crop(emi_voc_1,ire_extent)
emi_voc_2_ire<-crop(emi_voc_2,ire_extent)
emi_voc_3_ire<-crop(emi_voc_3,ire_extent)
emi_voc_4_ire<-crop(emi_voc_4,ire_extent)

coord_ire <- coordinates(emi_1_ire)

pred_results_nox_1_ire <- extract(emi_nox_1_ire,coord_ire)
pred_results_nox_2_ire <- extract(emi_nox_2_ire,coord_ire)
pred_results_nox_3_ire <- extract(emi_nox_3_ire,coord_ire)
pred_results_nox_4_ire <- extract(emi_nox_4_ire,coord_ire)
pred_results_voc_1_ire <- extract(emi_voc_1_ire,coord_ire)
pred_results_voc_2_ire <- extract(emi_voc_2_ire,coord_ire)
pred_results_voc_3_ire <- extract(emi_voc_3_ire,coord_ire)
pred_results_voc_4_ire <- extract(emi_voc_4_ire,coord_ire)


pred_results_nox_ire <- pred_results_nox_1_ire + pred_results_nox_2_ire + pred_results_nox_3_ire + pred_results_nox_4_ire
pred_results_voc_ire <- pred_results_voc_1_ire + pred_results_voc_2_ire + pred_results_voc_3_ire + pred_results_voc_4_ire
sum(is.na(pred_results_ire))
summary(pred_results_nox_ire)
summary(pred_results_voc_ire)

#reading output data
con_ire_somo <- con_somo[con_somo$CountryOrTerritory=="Ireland",]
sum(is.na(con_ire_somo))
summary(con_ire_somo)
max_ire_aq_somo<-max(con_ire_somo$AQValue)
#obtaining list of all cooridnates across Europe
coord_ire_1_somo <- cbind(con_ire_somo$SamplingPoint_Longitude, con_ire_somo$SamplingPoint_Latitude)

#extracting the data for required inputs 
results_1_nox <- extract(emi_nox_1_ire, coord_ire_1_somo)
results_2_nox <- extract(emi_nox_2_ire, coord_ire_1_somo)
results_3_nox <- extract(emi_nox_3_ire, coord_ire_1_somo)
results_4_nox <- extract(emi_nox_4_ire, coord_ire_1_somo)
results_1_voc <- extract(emi_voc_1_ire, coord_ire_1_somo)
results_2_voc <- extract(emi_voc_2_ire, coord_ire_1_somo)
results_3_voc <- extract(emi_voc_3_ire, coord_ire_1_somo)
results_4_voc <- extract(emi_voc_4_ire, coord_ire_1_somo)

#adding inputs from all results and cumulating from all macrosectors
final_results_ire_nox = results_1_nox+results_2_nox+results_3_nox+results_4_nox
final_results_ire_voc = results_1_voc+results_2_voc+results_3_voc+results_4_voc
summary(final_results_ire_voc)

#forming a new dataframe with coordinates,inputs and outputs
new_df_ire_somo <- data.frame(con_ire_somo$SamplingPoint_Longitude,con_ire_somo$SamplingPoint_Latitude,final_results_ire_nox,final_results_ire_voc,con_ire_somo$AQValue)
summary(new_df_ire_somo)
#converting data frame into spatial points data frame
spdf.con_ire.somo <- SpatialPointsDataFrame(coords = coord_ire_1_somo,data = new_df_ire_somo[,3:5],proj4string = crs(emi_2))
colnames(spdf.con_ire.somo@data)<-c("Emissions_nox","Emissions_voc","AQValue_somo")
spdf.pred_ire.somo <- SpatialPointsDataFrame(coords = coord_ire,data = as.data.frame(cbind(pred_results_nox_ire,pred_results_voc_ire)),proj4string = crs(emi_1))
colnames(spdf.pred_ire.somo@data)<-c("Emissions_nox","Emissions_voc")

#distance matrix
dmat_ire_somo<-gw.dist(dp.locat = coordinates(spdf.con_ire.somo))
#selecting bandwidth for geographically weighted regression
bw_optimal_bisq_ire_somo <- bw.gwr(AQValue_somo~Emissions_nox+Emissions_voc,data = spdf.con_ire.somo,approach = "CV",kernel = "bisquare",adaptive = TRUE,dMat = dmat_ire_somo)
#predicted results
gwr.pred_ire.somo<-gwr.predict(AQValue_somo~Emissions_nox+Emissions_voc,data=spdf.con_ire.somo,predictdata = spdf.pred_ire.somo,bw=bw_optimal_bisq_ire_somo,kernel = "bisquare",adaptive = TRUE, dMat1 = gw.dist(dp.locat = coordinates(spdf.con_ire.somo),rp.locat = coordinates(spdf.pred_ire.somo)), dMat2 = gw.dist(dp.locat = coordinates(spdf.con_ire.somo)))
pred_df_ire_somo<-as.data.frame(gwr.pred_ire.somo$SDF$prediction)
colnames(pred_df_ire_somo)<-"Pred_AQValue"
pred_df_ire_somo$Pred_AQValue[which(pred_df_ire_somo$Pred_AQValue<0)]=0
pred_df_ire_somo$Pred_AQValue[which(pred_df_ire_somo$Pred_AQValue>max_ire_aq_somo)]=max_ire_aq_somo
summary(pred_df_ire_somo)
hist(pred_df_ire_somo$Pred_AQValue)

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

f_plot<-ggplot(data=pred_df_ire_somo,aes(x=coord_ire[,1],y=coord_ire[,2]))
f_plot<-f_plot+geom_point(aes(colour=Pred_AQValue))+scale_color_gradient2(midpoint = 3000,low = "green",mid = "white",high = "red")+geom_polygon(data = counties,colour="black",fill=NA,aes(x=long, y=lat, group=group))
f_plot<-f_plot + labs(x="",y="",title = "Predicted Conc. (SOMO35) in Ireland")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
f_plot

#----------------------------------------------ITALY-------------------------------------------

it_extent<-extent(6.7499552751,18.4802470232, 36.619987291, 47.1153931748)
emi_nox_1_it<-crop(emi_nox_1,it_extent)
emi_nox_2_it<-crop(emi_nox_2,it_extent)
emi_nox_3_it<-crop(emi_nox_3,it_extent)
emi_nox_4_it<-crop(emi_nox_4,it_extent)
emi_voc_1_it<-crop(emi_voc_1,it_extent)
emi_voc_2_it<-crop(emi_voc_2,it_extent)
emi_voc_3_it<-crop(emi_voc_3,it_extent)
emi_voc_4_it<-crop(emi_voc_4,it_extent)

coord_it <- coordinates(emi_1_it)

pred_results_nox_1_it <- extract(emi_nox_1_it,coord_it)
pred_results_nox_2_it <- extract(emi_nox_2_it,coord_it)
pred_results_nox_3_it <- extract(emi_nox_3_it,coord_it)
pred_results_nox_4_it <- extract(emi_nox_4_it,coord_it)
pred_results_voc_1_it <- extract(emi_voc_1_it,coord_it)
pred_results_voc_2_it <- extract(emi_voc_2_it,coord_it)
pred_results_voc_3_it <- extract(emi_voc_3_it,coord_it)
pred_results_voc_4_it <- extract(emi_voc_4_it,coord_it)


pred_results_nox_it <- pred_results_nox_1_it + pred_results_nox_2_it + pred_results_nox_3_it + pred_results_nox_4_it
pred_results_voc_it <- pred_results_voc_1_it + pred_results_voc_2_it + pred_results_voc_3_it + pred_results_voc_4_it
sum(is.na(pred_results_it))
summary(pred_results_nox_it)
summary(pred_results_voc_it)

#reading output data
con_it_somo <- con_somo[con_somo$CountryOrTerritory=="Italy",]
sum(is.na(con_it_somo))
summary(con_it_somo)
max_it_aq_somo<-max(con_it_somo$AQValue)
#obtaining list of all cooridnates across Europe
coord_it_1_somo <- cbind(con_it_somo$SamplingPoint_Longitude, con_it_somo$SamplingPoint_Latitude)

#extracting the data for required inputs 
results_1_nox <- extract(emi_nox_1_it, coord_it_1_somo)
results_2_nox <- extract(emi_nox_2_it, coord_it_1_somo)
results_3_nox <- extract(emi_nox_3_it, coord_it_1_somo)
results_4_nox <- extract(emi_nox_4_it, coord_it_1_somo)
results_1_voc <- extract(emi_voc_1_it, coord_it_1_somo)
results_2_voc <- extract(emi_voc_2_it, coord_it_1_somo)
results_3_voc <- extract(emi_voc_3_it, coord_it_1_somo)
results_4_voc <- extract(emi_voc_4_it, coord_it_1_somo)

#adding inputs from all results and cumulating from all macrosectors
final_results_it_nox = results_1_nox+results_2_nox+results_3_nox+results_4_nox
final_results_it_voc = results_1_voc+results_2_voc+results_3_voc+results_4_voc
summary(final_results_it_nox)

#forming a new dataframe with coordinates,inputs and outputs
new_df_it_somo <- data.frame(con_it_somo$SamplingPoint_Longitude,con_it_somo$SamplingPoint_Latitude,final_results_it_nox,final_results_it_voc,con_it_somo$AQValue)
summary(new_df_it_somo)
#converting data frame into spatial points data frame
spdf.con_it.somo <- SpatialPointsDataFrame(coords = coord_it_1_somo,data = new_df_it_somo[,3:5],proj4string = crs(emi_2))
colnames(spdf.con_it.somo@data)<-c("Emissions_nox","Emissions_voc","AQValue_somo")
spdf.pred_it.somo <- SpatialPointsDataFrame(coords = coord_it,data = as.data.frame(cbind(pred_results_nox_it,pred_results_voc_it)),proj4string = crs(emi_1))
colnames(spdf.pred_it.somo@data)<-c("Emissions_nox","Emissions_voc")

#distance matrix
dmat_it_somo<-gw.dist(dp.locat = coordinates(spdf.con_it.somo))
#selecting bandwidth for geographically weighted regression
bw_optimal_bisq_it_somo <- bw.gwr(AQValue_somo~Emissions_nox+Emissions_voc,data = spdf.con_it.somo,approach = "CV",kernel = "bisquare",adaptive = TRUE,dMat = dmat_it_somo)
#predicted results
gwr.pred_it.somo<-gwr.predict(AQValue_somo~Emissions_nox+Emissions_voc,data=spdf.con_it.somo,predictdata = spdf.pred_it.somo,bw=bw_optimal_bisq_it_somo,kernel = "bisquare",adaptive = TRUE, dMat1 = gw.dist(dp.locat = coordinates(spdf.con_it.somo),rp.locat = coordinates(spdf.pred_it.somo)), dMat2 = gw.dist(dp.locat = coordinates(spdf.con_it.somo)))
pred_df_it_somo<-as.data.frame(gwr.pred_it.somo$SDF$prediction)
colnames(pred_df_it_somo)<-"Pred_AQValue"
pred_df_it_somo$Pred_AQValue[which(pred_df_it_somo$Pred_AQValue<0)]=0
pred_df_it_somo$Pred_AQValue[which(pred_df_it_somo$Pred_AQValue>max_it_aq_somo)]=max_it_aq_somo
summary(pred_df_it_somo)
hist(pred_df_it_somo$Pred_AQValue)

#plotting the predicted results for Ireland
spdf.shape_it<-readOGR("ITA_adm0.shp")
crs(spdf.shape_it)
spdf.shape_it@data$id<-rownames(spdf.shape_it@data)
spdf.points_it<-fortify(spdf.shape_it,region = "id")
states <- inner_join(spdf.points_it, spdf.shape_it@data, by="id")
#tm65 <- "+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=1.000035 +x_0=200000 +y_0=250000 +a=6377340.189 +b=6356034.447938534 +units=m +no_defs "
#newlonglat_ire<-project(cbind(counties$long,counties$lat),proj = tm65,inverse = TRUE)
#counties$long<-newlonglat_ire[,1]
#counties$lat<-newlonglat_ire[,2]

g_plot<-ggplot(data=pred_df_it_somo,aes(x=coord_it[,1],y=coord_it[,2]))
g_plot<-g_plot+geom_point(aes(colour=Pred_AQValue))+scale_color_gradient2(midpoint=6000,low = "green",mid = "white",high = "red")+geom_polygon(data = states,colour="black",fill=NA,aes(x=long, y=lat, group=group))
g_plot<-g_plot + labs(x="",y="",title = "Predicted Conc. (SOMO35) in Italy")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
g_plot
