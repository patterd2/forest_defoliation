#install.packages("MODISTools",dependencies=TRUE)
library(MODISTools)

mainDir="~/Rutgers/quantification_of_resilience/from_real_data/canada_spruce_budworm/location1/"

nyears=19

long_min<--68.23
long_max<--67.40

lat_min<-49.43
lat_max<-49.74
#############################################################################################################
#YOU CAN USE THIS PART TO FIND OUT WHETHER A SEPARATION OF 0.01 DD REPRESENTS MORE THAN 250m...##############
#install.packages("rgdal", dependencies = TRUE)
library(sp)
library(rgdal)

#Function
LongLatToUTM<-function(x,y,zone){
 xy <- data.frame(ID = 1:length(x), X = x, Y = y)
 coordinates(xy) <- c("X", "Y")
 proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
 res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
 return(as.data.frame(res))
}

zone<-"10V" #Valid for the narrow area we are considering.

#UTM_low_left_corner<-LongLatToUTM(long_min,lat_min,zone)
#UTM_low_right_corner<-LongLatToUTM(long_max,lat_min,zone)

#UTM_high_left_corner<-LongLatToUTM(long_min,lat_max,zone)
#UTM_high_right_corner<-LongLatToUTM(long_max,lat_max,zone)

#############################################################################################################


mt_products()

#According to Jepsen et al, cited by Vindstad et al, the data prouct they used is MODIS Terra NDVI with a 16-day temporal resolution, MOD13Q1. For each pixel, they focused on 4 16-day periods starting from day 177 within each year (i.e. late June to late Augst).

prod<-"MOD13Q1"

#day_init=177
#day_end=day_init+3*16   #It's 4 periods/points in total

mt_bands(prod)
colnames(mt_bands(prod))

band<-"250m_16_days_NDVI"

#File with the coordinates in decimal degrees (column 3 and 4):
#file_coordinates<-"../../source/Coordinates.Vindstad.txt"

folder_out<-"./"

#temp<-read.delim(file_coordinates)
#temp$name<-paste0(temp$Block,"_",temp$Transect)

longitudes<-seq(long_min,long_max,by=0.01)
latitudes<-seq(lat_min,lat_max,by=0.01)

temp<-expand.grid(longitudes,latitudes)
colnames(temp)<-c("X_dd","Y_dd")
temp$name<-paste0(match(temp$X_dd,longitudes),"_",match(temp$Y_dd,latitudes))

#Keep only X_dd and Y_dd, which are longitude and latitude, respectively, in decimal degrees
coordinates<-data.frame("site_name"=temp$name,"lat"=temp$Y_dd,"long"=temp$X_dd)  

ncoordinates<-nrow(coordinates)


##THIS IS JUST FOR ONE PIXEL. YOU CAN ALSO DEFINE LARGER REGIONS!!
##WE WILL DOWNLOAD ALL BANDS TO MAKE SURE THE PIXEL QUALITY/RELIABILITY DATA ARE THERE AS WELL
for (j in 1:nyears)
{
number<-sprintf("%02.0f",j)

year<-1999+as.numeric(number)
print(year)

subDir=as.character(year)

dir.create(file.path(mainDir, subDir))
setwd(file.path(mainDir, subDir))

#date_first<-as.Date(day_init, origin = paste0(year,"-01-01"))
#date_last<-as.Date(day_end, origin = paste0(year,"-01-01"))
date_first<-paste0(year,"-01-01")
date_last<-paste0(year,"-12-31")

#mt_batch_subset(df = coordinates, product = prod, band = band, start = date_first, end = date_last, km_lr = 0, km_ab = 0, out_dir = folder_out, internal = FALSE)

#name<-paste0("subset_",year)
#assign(name,mt_batch_subset(df = coordinates, product = prod, band = band, start = date_first, end = date_last, km_lr = 0, km_ab = 0, internal = TRUE))
#mt_write(df = name,out_dir = folder_out)

#If you have problems with the batch function saturating connection sockets, try instead:
for (k in 1:ncoordinates){

destfile<-paste0(coordinates[k,1],"_",prod,"_",date_first,"_",date_last,".csv")
if(!file.exists(destfile)){
mt_subset(product = prod, lat=coordinates[k,2], lon=coordinates[k,3], site_name=as.character(coordinates[k,1]), band = NULL, start = date_first, end = date_last, km_lr = 0, km_ab = 0, out_dir = folder_out, internal = FALSE)
}
info = file.info(destfile)
if(info$size == 0){
mt_subset(product = prod, lat=coordinates[k,2], lon=coordinates[k,3], site_name=as.character(coordinates[k,1]), band = NULL, start = date_first, end = date_last, km_lr = 0, km_ab = 0, out_dir = folder_out, internal = FALSE)
}
}
}

#Test/example: first pixel on the coordinates list, last year
#subset<-mt_subset(product = prod,lat = 70.03051605,lon = 28.55617521,band = band,start = date_first,end = date_last, internal = TRUE)


##TEMPERATURE DATA FOR THE SAME PIXELS/LOCATIONS:

prod<-"MOD11A2"
mt_bands(prod)  

band<-"LST_Day_1km"

for (j in 1:nyears)
{
number<-sprintf("%02.0f",j)

year<-1999+as.numeric(number)
print(year)

subDir=as.character(year)

#dir.create(file.path(mainDir, subDir))
setwd(file.path(mainDir, subDir))

#date_first<-as.Date(day_init, origin = paste0(year,"-01-01"))
#date_last<-as.Date(day_end, origin = paste0(year,"-01-01"))
date_first<-paste0(year,"-01-01")
date_last<-paste0(year,"-12-31")

for (k in 1:ncoordinates){
destfile<-paste0(coordinates[k,1],"_",prod,"_",date_first,"_",date_last,".csv")
if(!file.exists(destfile)){
mt_subset(product = prod, lat=coordinates[k,2], lon=coordinates[k,3], site_name=as.character(coordinates[k,1]), band = band, start = date_first, end = date_last, km_lr = 0, km_ab = 0, out_dir = folder_out, internal = FALSE)
}
info = file.info(destfile)
if(info$size == 0){
mt_subset(product = prod, lat=coordinates[k,2], lon=coordinates[k,3], site_name=as.character(coordinates[k,1]), band = NULL, start = date_first, end = date_last, km_lr = 0, km_ab = 0, out_dir = folder_out, internal = FALSE)
}
}

}
