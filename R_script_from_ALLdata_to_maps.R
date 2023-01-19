library(zoo)
library(bfast)
library(lubridate)
library(tibble) #To add columns to data frames
library(crone) #For the heaviside function
library(ks)
#library(logspline)

library("LaplacesDemon") #To use the Modes function
library("fANCOVA") #For the optimized loess parametrization function

library(rootSolve)

library(gplots)

library("ggplot2")


mytheme <-
  theme(
    plot.title = element_text(
      family = "Helvetica",
      face = "bold",
      size = (15)
    ),
    legend.title = element_text(
      colour = "steelblue",
      face = "bold.italic",
      family = "Helvetica"
    ),
    legend.text = element_text(face = "italic", family = "Helvetica"),
    axis.title = element_text(family = "Helvetica", size = (16)),
    axis.text = element_text(family = "Helvetica", size = (12)),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(
      size = 0.5,
      linetype = "solid",
      colour = "black"
    )
    #                  panel.background = element_rect(fill = "white",colour = "white",size = 0.0, linetype = "solid"),
    #                  panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "green"),
    #                  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "red")
  )

#For the map part:
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

# load library
library(ggmap)
###

step_f <- function(x) {
  step_f <- x * heaviside(x, x0 = 0.)
  if (is.na(x) == TRUE) {
    step_f <- 0
  }
  return(step_f)
}

step_f2 <- function(x) {
  step_f2 <- x * heaviside(x, x0 = 0.)
  if (is.na(x) == TRUE) {
    step_f2 <- 0
  } else{
    if (x > 1) {
      step_f2 <- 1
    }
  }
  return(step_f2)
}

step_f3 <- function(x) {
  step_f3 <- x * heaviside(x, x0 = 0.)
  return(step_f3)
}

mode_f <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

first_year <-
  2001#Because 2000 was a weird year for dates when data were taken that screws up the periodicity.
year0 <- first_year - 1
last_year <- 2018

#Period for the stable/normal behavior:
normal_start <- last_year - 1
normal_end <- last_year + 1

nyears = last_year - year0

setwd(
  "/Users/Denis/Downloads/forest_defoliation/canada_spruce_budworm/location1"
)
mainDir = "./"

NDVI_threshold <-
  0.25#Threshold under which we consider that's not even a tree (see when used below for more details).

summer_init <- 177
summer_end <- summer_init + 3 * 16   #It's 4 periods/points in total

summerM_init <- 5
summerM_end <- 10

Tstart <- 1
Tend <- summer_init

Tlow <-
  5 + 273.15 #Typical lower T threshold needed for development to occur, in Kelvin.
Thigh <-
  40 + 273.15 #We don't have info for most cases, so we just make it really high.

long_min <- -68.23
long_max <- -67.40

lat_min <- 49.43
lat_max <- 49.74

longitudes <- seq(long_min, long_max, by = 0.01)
latitudes <- seq(lat_min, lat_max, by = 0.01)

temp <- expand.grid(longitudes, latitudes)
colnames(temp) <- c("X_dd", "Y_dd")
temp$name <-
  paste0(match(temp$X_dd, longitudes), "_", match(temp$Y_dd, latitudes))

#Keep only X_dd and Y_dd, which are longitude and latitude, respectively, in decimal degrees
coordinates <-
  data.frame(
    "site_name" = temp$name,
    "lat" = temp$Y_dd,
    "long" = temp$X_dd
  )

ncoordinates <- nrow(coordinates)

prod_T <- "MOD11A2"
prod_NDVI <- "MOD13Q1"
prod_TC <- "MOD44B"
prod_LU <- "MCD12Q1"

band_T <-
  "LST_Day_1km"#No need because it's the only band we downloaded!!
band_rel <- "250m_16_days_pixel_reliability"
band_QA <- "250m_16_days_VI_Quality"
band_NDVI <- "250m_16_days_NDVI"

band_TC <- "Percent_Tree_Cover"
band_NOTC <- "Percent_NonTree_Vegetation"
band_NOVEG <- "Percent_NonVegetated"
band_QATC <- "Quality"

band_LUT1 <- "LC_Type1"
band_LUT2 <- "LC_Type2"
band_LUT3 <- "LC_Type3"
band_LUT4 <- "LC_Type4"
band_LUT5 <- "LC_Type5"
band_LUP1 <- "LC_Prop1"
band_LUP2 <- "LC_Prop2"
band_LUP3 <- "LC_Prop3"
band_LUP1conf <- "LC_Prop1_Assessment"
band_LUP2conf <- "LC_Prop2_Assessment"
band_LUP3conf <- "LC_Prop3_Assessment"
band_LUQA <- "QC"
band_LULW <- "LW"

#According to user's manual and the column "scale_factor", NDVI has a scale factor 1/10000, i.e. you need to multiply the stored number by that factor to get the real NDVI score. See https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13q1
scale <- 0.0001
scaleT <- 0.02#In Kelvin, check out any of the csv files!!

folder_out <- paste0("./final_files_DEF_FOR_REAL/")
dir.create(file.path(mainDir, folder_out))

DEF_ALL <-
  data.frame(
    "site_name" = character(),
    "lat_dd" = character(),
    "long_dd" = character(),
    "year" = integer(),
    "date" = as.Date(character()),
    "day" = integer(),
    "temp_Kelvin" = double(),
    "NDVI" = double(),
    "pixel_rel" = integer(),
    "QA" = integer(),
    dec_date_ADJ = double(),
    prediction = double(),
    dec_date_NDVImax = double(),
    NDVImax = double(),
    rel_diff = double()
  )

DEF_ALL_w_NAs <-
  data.frame(
    "site_name" = character(),
    "lat_dd" = character(),
    "long_dd" = character(),
    "year" = integer(),
    "date" = as.Date(character()),
    "day" = integer(),
    "temp_Kelvin" = double(),
    "NDVI" = double(),
    "pixel_rel" = integer(),
    "QA" = integer(),
    dec_date_ADJ = double(),
    prediction = double(),
    dec_date_NDVImax = double(),
    NDVImax = double(),
    rel_diff = double()
  )

for (i in 2469:ncoordinates) {
  #for every location...
  
  site_name <- as.character(coordinates[i, 1])
  
  print(paste(i, "out of", ncoordinates, "=", site_name))
  
  data <-
    data.frame(
      "site_name" = character(),
      "lat_dd" = character(),
      "long_dd" = character(),
      "year" = integer(),
      "date" = as.Date(character()),
      "day" = integer(),
      "temp_Kelvin" = double(),
      "NDVI" = double(),
      "pixel_rel" = integer(),
      "QA" = integer()
    )
  
  dataT <-
    data.frame(
      "site_name" = character(),
      "lat_dd" = character(),
      "long_dd" = character(),
      "year" = integer(),
      "date" = as.Date(character()),
      "day" = integer(),
      "temp_Kelvin" = double()
    )
  
  for (k in 1:nyears) {
    #For every year...
    
    number <- sprintf("%02.0f", k)
    year <- year0 + as.numeric(number)
    #print(year)
    date_first <- paste0(year, "-01-01")
    date_last <- paste0(year, "-12-31")
    
    folder_data <- paste0("./data/", year, "/")
    
    file_inT <-
      paste0(folder_data,
             site_name,
             "_",
             prod_T,
             "_",
             date_first,
             "_",
             date_last,
             ".csv")
    file_inN <-
      paste0(folder_data,
             site_name,
             "_",
             prod_NDVI,
             "_",
             date_first,
             "_",
             date_last,
             ".csv")
    file_inTC <-
      paste0(folder_data,
             site_name,
             "_",
             prod_TC,
             "_",
             date_first,
             "_",
             date_last,
             ".csv")
    file_inLU <-
      paste0(folder_data,
             site_name,
             "_",
             prod_LU,
             "_",
             date_first,
             "_",
             date_last,
             ".csv")
    
    #Read the file into a dataframe, add it to a temporary dataframe, ignoring the comments (as most of that info is in the name file)
    
    subsetN <-
      read.csv(
        file_inN,
        header = TRUE,
        fill = TRUE,
        sep = ",",
        comment.char = '#'
      ) # ERROR HERE: DP - Jan 19
    
    #If the entry for the date is not in the temperature, tree cover, or land cover files, add the entry to their data frames using a NA value:
    
    if (!file.exists(file_inT)) {
      #subsetT<-data.frame(band=c(as.character(band_T)),data=NA)
      subsetT <- subsetN[subsetN$band == band_NDVI, ]
      subsetT$band <- band_T
      subsetT$data <- NA
    } else{
      subsetT <-
        read.csv(
          file_inT,
          header = TRUE,
          fill = TRUE,
          sep = ",",
          comment.char = '#'
        )
    }
    
    if (!file.exists(file_inTC)) {
      subsetTC <-
        data.frame(
          band = c(
            as.character(band_TC),
            as.character(band_NOTC),
            as.character(band_NOVEG),
            as.character(band_QATC)
          ),
          data = NA
        )
    } else{
      subsetTC <-
        read.csv(
          file_inTC,
          header = TRUE,
          fill = TRUE,
          sep = ",",
          comment.char = '#'
        )
    }
    
    if (!file.exists(file_inLU)) {
      subsetLU <-
        data.frame(
          band = c(
            as.character(band_LUT1),
            as.character(band_LUT2),
            as.character(band_LUT3),
            as.character(band_LUT4),
            as.character(band_LUT5),
            as.character(band_LUP1),
            as.character(band_LUP2),
            as.character(band_LUP3),
            as.character(band_LUP1conf),
            as.character(band_LUP2conf),
            as.character(band_LUP3conf),
            as.character(band_LUQA),
            as.character(band_LULW)
          ),
          data = NA
        )
    } else{
      subsetLU <-
        read.csv(
          file_inLU,
          header = TRUE,
          fill = TRUE,
          sep = ",",
          comment.char = '#'
        )
    }
    
    #Keep only the bands we are interested in:
    
    tempT <- subsetT[subsetT$band == band_T, ]
    tempNDVIpre <- subsetN[subsetN$band == band_NDVI, ]
    temprel <- subsetN[subsetN$band == band_rel, ]
    tempQA <- subsetN[subsetN$band == band_QA, ]
    
    tempTC <- subsetTC[subsetTC$band == band_TC, ]
    tempNOTC <- subsetTC[subsetTC$band == band_NOTC, ]
    tempNOVEG <- subsetTC[subsetTC$band == band_NOVEG, ]
    tempQATC <- subsetTC[subsetTC$band == band_QATC, ]
    
    tempLU <- subsetLU
    
    tempT$data <- as.numeric(as.character(tempT$data))
    tempNDVIpre$data <- as.numeric(as.character(tempNDVIpre$data))
    temprel$data <- as.numeric(as.character(temprel$data))
    
    tempQA$data <- as.numeric(as.character(tempQA$data))
    tempTC$data <- as.numeric(as.character(tempTC$data))
    tempNOTC$data <- as.numeric(as.character(tempNOTC$data))
    tempNOVEG$data <- as.numeric(as.character(tempNOVEG$data))
    tempQATC$data <- as.numeric(as.character(tempQATC$data))
    
    tempLU$data <- as.numeric(as.character(tempLU$data))
    
    tempNDVI <-
      data.frame(
        "day" = yday(tempNDVIpre$calendar_date),
        tempNDVIpre[c("modis_date",
                      "calendar_date",
                      "tile",
                      "proc_date",
                      "pixel")],
        "NDVI" = tempNDVIpre$data,
        "pixel_rel" = temprel$data,
        "QA" = tempQA$data
      )
    
    tempdata <-
      data.frame(
        "site_name" = site_name,
        "lat_dd" = coordinates[i, 2],
        "long_dd" = coordinates[i, 3],
        "year" = year,
        "date" = as.Date(tempNDVI$calendar_date),
        "day" = yday(tempNDVI$calendar_date),
        "NDVI" = scale * tempNDVI$NDVI,
        "pixel_rel" = tempNDVI$pixel_rel,
        "QA" = tempQA$data
      )
    data <- rbind(data, tempdata)
    #file_out<-paste0(folder_out,"NDVI_and_T-",year,"-site_",site_name,".txt")
    #write.table(temp_df,file_out, sep="\t",  row.names=FALSE, col.names=TRUE)
    
    tempdata <-
      data.frame(
        "site_name" = site_name,
        "lat_dd" = coordinates[i, 2],
        "long_dd" = coordinates[i, 3],
        "year" = year,
        "date" = as.Date(tempT$calendar_date),
        "day" = yday(tempT$calendar_date),
        "temp_Kelvin" = scaleT * tempT$data
      )
    dataT <- rbind(dataT, tempdata)
    
  }#End of the years loop, i.e. we have all the data for this location across all the years.
  
  #Early test to spot water locations: if most NDVI or T data are invalid:
  
  nBADT <- length(which(dataT$temp_Kelvin == 0.))
  nBADNDVI <- length(which(data$pixel_rel == -1))
  if ((nBADT < 0.25 * nrow(dataT) && nBADNDVI < 0.25 * nrow(data))) {
    #Replace all T=0s, which are artifacts, with NaNs (which is different from NA, which are included because of missing data or dates in which NDVI was measured but T was not, and which is an important distinction for the time series and model fit):
    
    dataT$temp_Kelvin[which(dataT$temp_Kelvin == 0.)] <- NaN
    
    
    #SELECTION OF "STANDARD NDVI" BEHAVIOR/MODEL##########################################
    
    NDVIts <- bfastts(data$NDVI, data$date, type = "16-day")
    
    #By default, we use the year with the second-highest SUMMER NDVI of all as representative of "normal" behavior (as in the original moth paper).
    #First, choose the highest NDVI of each summer:
    tempsummer <-
      data[which(data$day %in% seq(summer_init, summer_end, 1)), ]
    y = first_year
    j = 1
    tempMAX <- 0
    tempy <- 0
    while (y <= last_year) {
      temp <- tempsummer[which(tempsummer$year == y), ]
      tempMAX[j] <- max(temp$NDVI)
      tempy[j] <- y
      y = y + 1
      j = j + 1
    }
    
    #Then choose the second best summer:
    jMAX2 <- order(tempMAX, decreasing = TRUE)[2]
    MAX2year <- tempy[jMAX2]
    
    #As long as the data for the whole year are as "complete" as possible so that you can do a proper fitting:
    temp_loc <- which(data$year == MAX2year)
    ntemploc <- length(temp_loc)
    j = 3
    while (ntemploc < 20 & j < nyears) {
      jMAX2 <- order(tempMAX, decreasing = TRUE)[j]
      MAX2year <- tempy[jMAX2]
      #As long as the data for the whole year are as "complete" as possible so that you can do a proper fitting:
      temp_loc <- which(data$year == MAX2year)
      ntemploc <- length(temp_loc)
      j + 1
    }
    
    yhistory = MAX2year
    ystart = yhistory + 1
    
    #YOU SHOULD CHECK WHETHER CONCLUSIONS CHANGE IF YOU TAKE THE YEAR OF THE MAX:
    ##Then choose the best summer:
    #jMAX<-order(tempMAX,decreasing=TRUE)[1]
    #MAXyear<-tempy[jMAX]
    #
    ##As long as the data for the whole year are as "complete" as possible so that you can do a proper fitting:
    #temp_loc<-which(data$year==MAXyear)
    #ntemploc<-length(temp_loc)
    #j=3
    #while (ntemploc<20 & j<nyears){
    #jMAX<-order(tempMAX,decreasing=TRUE)[j]
    #MAXyear<-tempy[jMAX]
    ##As long as the data for the whole year are as "complete" as possible so that you can do a proper fitting:
    #temp_loc<-which(data$year==MAXyear)
    #ntemploc<-length(temp_loc)
    #j+1
    #}
    #
    #OR THAT OF THE CLOSEST TO THE MEDIAN (BUT IMMEDIATELY ABOVE).
    #tempsummer<-data[which(data$day %in% seq(summer_init,summer_end,1)),]
    #dist_to_median<-1.-tempMAX/median(tempMAX)
    #dist_to_medianNEG<-dist_to_median[which(dist_to_median<=0)]
    #j=order(dist_to_medianNEG,decreasing=TRUE)[1]
    #dist<-dist_to_medianNEG[j]
    #yhistory<-tempy[which(dist_to_median==dist)]
    #ystart=yhistory+1
    
    #Fit a harmonic model to the stable period:
    #monshortalt <- bfastmonitor(NDVIts, start = ystart,history=yhistory)
    monshortalt <-
      bfastmonitor(
        NDVIts,
        formula = response ~ harmon,
        start = ystart ,
        history = yhistory
      )
    #plot(monshortalt)
    
    if (length(monshortalt$model$coefficients) == 8) {
      #print("We are using the trend too")
      
      freq <- frequency(NDVIts)
      f <- function(x) {
        #j<-which(names(monshortalt$model$residuals)==as.character(x))
        imax <- max(as.numeric(rownames(monshortalt$model$model)))
        imin <- min(as.numeric(rownames(monshortalt$model$model)))
        #j<-which(seq(imin,imax)==x)
        #y=j-1.0
        y = (x - imin)
        while (y < 0) {
          y = y + freq
        }
        f <-
          monshortalt$model$coefficients[1][[1]] + x * monshortalt$model$coefficients[2][[1]] +
          monshortalt$model$coefficients[3][[1]] * cos(2 * pi * 1 * y / freq) + monshortalt$model$coefficients[4][[1]] *
          cos(2 * pi * 2 * y / freq) + monshortalt$model$coefficients[5][[1]] * cos(2 *
                                                                                      pi * 3 * y / freq) + monshortalt$model$coefficients[6][[1]] * sin(2 * pi *
                                                                                                                                                          1 * y / freq) + monshortalt$model$coefficients[7][[1]] * sin(2 * pi * 2 *
                                                                                                                                                                                                                         y / freq) + monshortalt$model$coefficients[8][[1]] * sin(2 * pi * 3 * y /
                                                                                                                                                                                                                                                                                    freq)
        print("THIS ONE IS NOT GOING TO WOOOOORK!!!!")
        return(f)
      }
      
    } else{
      #print("No trend considered")
      
      freq <- frequency(NDVIts)
      f <- function(x) {
        #j<-which(names(monshortalt$model$residuals)==as.character(x))
        imax <- max(as.numeric(rownames(monshortalt$model$model)))
        imin <- min(as.numeric(rownames(monshortalt$model$model)))
        #j<-which(seq(imin,imax)==x)
        #y=j-1.0
        y = (x - imin)
        while (y < 0) {
          y = y + freq
        }
        f <-
          monshortalt$model$coefficients[1][[1]] + monshortalt$model$coefficients[2][[1]] *
          cos(2 * pi * 1 * y / freq) + monshortalt$model$coefficients[3][[1]] * cos(2 *
                                                                                      pi * 2 * y / freq) + monshortalt$model$coefficients[4][[1]] * cos(2 * pi *
                                                                                                                                                          3 * y / freq) + monshortalt$model$coefficients[5][[1]] * sin(2 * pi * 1 *
                                                                                                                                                                                                                         y / freq) + monshortalt$model$coefficients[6][[1]] * sin(2 * pi * 2 * y /
                                                                                                                                                                                                                                                                                    freq) + monshortalt$model$coefficients[7][[1]] * sin(2 * pi * 3 * y / freq)
        return(f)
      }
    }
    
    temp <- seq(1, length(NDVIts))
    ilist <- temp[!(temp %in% which(is.na(NDVIts)))]
    model <- unlist(lapply(ilist, f))
    
    
    data$dec_date_ADJ <- time(NDVIts)[which(!is.na(NDVIts))]
    data$prediction <- model
    data$year_ref <- yhistory
    
    
    file <- paste0(folder_out, "NDVIseries_site_", site_name, ".png")
    png(file)
    plot(
      data$dec_date_ADJ,
      data$NDVI,
      xlab = "Year",
      ylab = "NDVI",
      main = "",
      pch = 20,
      type = "l",
      lwd = 2,
      lty = 1,
      col = "blue",
      xlim = c(first_year, last_year),
      xaxt = 'n'
    )
    axis(1, at = seq(2000, 2019, by = 0.5), las = 2)
    par(new = TRUE)
    lines(
      data$dec_date_ADJ,
      data$prediction,
      lwd = 2,
      lty = 1,
      col = "red",
      xlim = c(first_year, last_year)
    )
    dev.off()
    
    ###########################################################################################
    
    #############CONSTRUCTION OF THE 8-day-FREQUENCY TEMPERATURE DATA FRAME AND ANNUAL REPRESENTATIVE################
    
    #Add any missing date that breaks the 8-day periodicity/frequency of data collection:
    year = first_year
    while (year <= last_year) {
      temp1 <-
        seq.Date(as.Date(paste0(year, "-01-01")), as.Date(paste0(year, "-12-31")), 8)
      temp2 <- dataT$date[which(dataT$year == year)]
      missing <- temp1[!(temp1 %in% temp2)]
      j = 1
      while (j <= length(missing)) {
        #print(missing[j])
        temp <- dataT[1, ]
        temp$date <- missing[j]
        temp$day <- yday(temp$date)
        temp$temp_Kelvin <- NA
        dataT <- rbind(dataT, temp)
        j = j + 1
      }
      year = year + 1
    }
    dataT <- dataT[order(dataT$date), ]
    rownames(dataT) <- seq(length = nrow(dataT))
    
    #Add the decimal version of the date:
    freq <- 2 * frequency(NDVIts)
    j = 0
    year = first_year
    temp <- 0
    while (year <= last_year) {
      temp[seq(j + 1, j + freq)] <-
        seq(year, year + 1 - 1 / freq, by = 1 / freq)
      j = j + freq
      year = year + 1
    }
    dataT$dec_date_ADJ <- temp
    
    #CONSTRUCT A TIME SERIES WITH dataT (NOW THAT IT'S BEEN BASICALLY MADE ONE):
    temp <- dataT
    #temp$temp_Kelvin[is.nan(temp$temp_Kelvin)]<-NA
    Tts <- ts(temp$temp_Kelvin,
              start = c(first_year, 1),
              freq = 46)
    
    
    #REFERENCE: YEAR WE CHOOSE AS REFERENCE FOR HEALTHY NDVI, UNDER THE ASSUMPTION THAT ALL CONDITIONS (TEMPERATURE AND INSECT BEHAVIOR) ARE ALSO IDEAL FOR NDVI THAT YEAR:
    
    
    yhistory = data$year_ref[1]
    ystart = yhistory + 1
    
    #Fit a harmonic model to the stable period:
    #monshortalt <- bfastmonitor(Tts, start = ystart,history=yhistory)
    monshortalt <-
      bfastmonitor(
        Tts,
        formula = response ~ harmon,
        start = ystart ,
        history = yhistory
      )
    #plot(monshortalt)
    
    if (length(monshortalt$model$coefficients) == 8) {
      #print("We are using the trend too")
      
      freq <- frequency(Tts)
      f <- function(x) {
        #j<-which(names(monshortalt$model$residuals)==as.character(x))
        imax <- max(as.numeric(rownames(monshortalt$model$model)))
        imin <- min(as.numeric(rownames(monshortalt$model$model)))
        #j<-which(seq(imin,imax)==x)
        #y=j-1.0
        y = (x - imin)
        while (y < 0) {
          y = y + freq
        }
        f <-
          monshortalt$model$coefficients[1][[1]] + x * monshortalt$model$coefficients[2][[1]] +
          monshortalt$model$coefficients[3][[1]] * cos(2 * pi * 1 * y / freq) + monshortalt$model$coefficients[4][[1]] *
          cos(2 * pi * 2 * y / freq) + monshortalt$model$coefficients[5][[1]] * cos(2 *
                                                                                      pi * 3 * y / freq) + monshortalt$model$coefficients[6][[1]] * sin(2 * pi *
                                                                                                                                                          1 * y / freq) + monshortalt$model$coefficients[7][[1]] * sin(2 * pi * 2 *
                                                                                                                                                                                                                         y / freq) + monshortalt$model$coefficients[8][[1]] * sin(2 * pi * 3 * y /
                                                                                                                                                                                                                                                                                    freq)
        print("THIS ONE IS NOT GOING TO WORK!!!!")
        return(f)
      }
      
    } else{
      #print("No trend considered")
      
      freq <- frequency(Tts)
      f <- function(x) {
        #j<-which(names(monshortalt$model$residuals)==as.character(x))
        imax <- max(as.numeric(rownames(monshortalt$model$model)))
        imin <- min(as.numeric(rownames(monshortalt$model$model)))
        #j<-which(seq(imin,imax)==x)
        #y=j-1.0
        y = (x - imin)
        while (y < 0) {
          y = y + freq
        }
        f <-
          monshortalt$model$coefficients[1][[1]] + monshortalt$model$coefficients[2][[1]] *
          cos(2 * pi * 1 * y / freq) + monshortalt$model$coefficients[3][[1]] * cos(2 *
                                                                                      pi * 2 * y / freq) + monshortalt$model$coefficients[4][[1]] * cos(2 * pi *
                                                                                                                                                          3 * y / freq) + monshortalt$model$coefficients[5][[1]] * sin(2 * pi * 1 *
                                                                                                                                                                                                                         y / freq) + monshortalt$model$coefficients[6][[1]] * sin(2 * pi * 2 * y /
                                                                                                                                                                                                                                                                                    freq) + monshortalt$model$coefficients[7][[1]] * sin(2 * pi * 3 * y / freq)
        return(f)
      }
    }
    
    #ADDING THE MODEL TO dataT, REPLACING WITH PREDICTIONS THE EXISTING NaNs or NAs:
    m = 1
    modelT <- NA
    while (m <= length(Tts)) {
      modelT[m] <- f(m)
      m = m + 1
    }
    
    dataT$predictionT <- modelT
    dataT$year_refT <- yhistory
    
    #ADDING THE MODEL TO THE NDVI DATA FRAME; SOME ADJUSTMENTS ARE NEEDED DUE TO THE DIFFERENCE IN FREQUENCY: given that both data and dataT have the "regularized" dec_date_ADJs, and that the latter's frequency is exactly twice the former, the former elements will always be included in the latter's list, which makes it easy to just pick from dataT the right elements:
    
    m = 1
    modelT <- 0
    realT <- 0
    while (m <= length(data$dec_date_ADJ)) {
      loc <- which.min(abs(data$dec_date_ADJ[m] - dataT$dec_date_ADJ))
      modelT[m] <- dataT$predictionT[loc]
      realT[m] <- dataT$temp_Kelvin[loc]
      m = m + 1
    }
    
    data$temp_Kelvin <- realT
    data$predictionT <- modelT
    data$year_refT <- yhistory
    
    
    #THE FOLLOWING PLOT SHOULD SHOW HOW NON-PERIODIC TEMPERATURE IS (SOME YEARS WARM UP SOONER THAN OTHERS, ETC!):
    file <- paste0(folder_out, "Tseries_site_", site_name, ".png")
    png(file)
    plot(
      dataT$dec_date_ADJ,
      dataT$temp_Kelvin,
      xlab = "Year",
      ylab = "Temperature (Kelvin)",
      main = "",
      pch = 20,
      type = "l",
      lwd = 2,
      lty = 1,
      col = "blue",
      xlim = c(first_year, last_year),
      xaxt = 'n'
    )
    axis(1, at = seq(2000, 2019, by = 0.5), las = 2)
    par(new = TRUE)
    lines(
      dataT$dec_date_ADJ,
      dataT$predictionT,
      lwd = 2,
      lty = 1,
      col = "red",
      xlim = c(first_year, last_year)
    )
    dev.off()
    
    
    #########################################################
    #YEARLY REPRESENTATIVE FOR T: THE DEGREE DAY ACCUMULATED UNTIL SUMMER, USING AS A REFERENCE THE YEAR WITH THE (SECOND) HIGHEST NDVI (AS INDICATOR OF NO INSECT ACTIVITY):
    
    #FOR EACH YEAR...
    
    j = 1
    year = first_year
    degree_day <- 0
    while (year <= last_year) {
      #print(year)
      
      #1) FIT A CONTINUOUS CURVE FOR TEMPERATURE
      temp_loc <- which(!is.na(dataT$temp_Kelvin) & dataT$year == year)
      tempT <- dataT$temp_Kelvin[temp_loc]
      temp_date <- dataT$dec_date_ADJ[temp_loc]
      
      fit <-
        smooth.spline(temp_date, tempT, cv = FALSE)#Find a smooth curve resembling the points
      fit_func <- splinefun(temp_date, fit$y)#Make a function out of it
      
      ###############
      ######ALTERNATIVE THAT IS QUALITATIVELY SIMILAR (BUT IS USED MORE OFTEN):
      #2) SWEEP EVERY DAY AND ADD THE NUMBER OF DEGREES-DAY ABOVE Tlow (IF YOU KNOW Thigh, MAYBE YOU SHOULD ALSO ADD THE CONDITION T<Thigh FOR ACCUMULATION):
      #days<-seq(year,year+1,1./365)#You can refine more, checking how many days exactly there was that year (if it's a leap year).
      days <- seq(year, year + summer_init / 365, 1. / 365)#Only until the summer
      
      degree_day[j] <- sum(step_f3(fit_func(days) - Tlow))
      
      j = j + 1
      year = year + 1
    }
    #plot(temp_date,fit_func(temp_date))
    #lines(days,fit_func(days))
    
    
    Tyear_ref <- unique(data$year_ref) #This should be only one year.
    
    #########################################################
    
    
    #In days within the year:
    summer_only <-
      data[which(data$day %in% seq(summer_init, summer_end, 1)), ]
    
    #FILTERING DATA FOR QUALITY:
    summer_only <-
      summer_only[summer_only$pixel_rel >= 0, ]#To filter the -1, which are fill-in values (also characterized by a NDVI of -0.3)
    #summer_only<-summer_only[summer_only$pixel_rel==0,]
    summer_only <-
      summer_only[summer_only$pixel_rel < 2, ] #This one is less strict than ==0, because it allows 0s and 1s; if we use this one, we have to remove artifacts like the sudden dips in NDVI that happen in one single measurement of the summer.
    summer_only <-
      summer_only[which(!(
        summer_only$NDVI < NDVI_threshold &
          summer_only$pixel_rel > 0
      )), ]#One way to remove the biggest artifacts.
    
    #Choose a representative per year:
    j = 1
    yearly_rep <-
      data.frame(matrix(ncol = length(names(summer_only)), nrow = 0))
    colnames(yearly_rep) <- names(summer_only)
    yearly_rep <-
      add_column(yearly_rep,
                 dec_date_NDVImax = NA,
                 .after = "year_refT")
    yearly_rep <-
      add_column(yearly_rep, NDVImax = NA, .after = "dec_date_NDVImax")
    yearly_rep <- add_column(yearly_rep, degree_day = NA, .after = "NDVImax")
    yearly_rep <- add_column(yearly_rep, Tyear_ref = NA, .after = "degree_day")
    
    year = first_year
    while (year <= last_year) {
      temp_loc <- which(summer_only$year == year)
      ntemploc <- length(temp_loc)
      if (ntemploc > 0) {
        temp <- summer_only$prediction[temp_loc]
        loc_w_max <- temp_loc[which(temp == max(temp))]
        yearly_rep[j, ] <- summer_only[loc_w_max, ]
        
        temp_loc <- which(floor(summer_only$dec_date_ADJ) == year)
        temp <- summer_only$NDVI[temp_loc]
        loc_w_max <- temp_loc[which(temp == max(temp))]
        if (length(loc_w_max) > 1) {
          loc_w_max <-
            loc_w_max[1]
        }#This happens if more than one location shows the max, in which case it doesn't matter which one we choose.
        #If that NDVImax<NDVI_threshold, the site is most probably not covered by trees:
        if (summer_only$NDVI[loc_w_max] < NDVI_threshold) {
          print(paste(
            year,
            ": no reliable representative for the year because no trees?!!"
          ))
          yearly_rep$dec_date_NDVImax[j] <- NA
          yearly_rep$NDVImax[j] <- NA
          yearly_rep$degree_day[j] <- NA
        } else{
          yearly_rep$dec_date_NDVImax[j] <- summer_only$dec_date_ADJ[loc_w_max]
          yearly_rep$NDVImax[j] <- summer_only$NDVI[loc_w_max]
        }
      } else{
        print(paste(year, ": no reliable representative for the year!!"))
        temp_loc <- which(data$year == year)
        temp <- data$prediction[temp_loc]
        loc_w_max <- temp_loc[which(temp == max(temp))]
        yearly_rep[j, ] <- data[loc_w_max, ]
        
        yearly_rep$dec_date_NDVImax[j] <- NA
        yearly_rep$NDVImax[j] <- NA
        yearly_rep$degree_day[j] <- NA
        
      }
      
      j = j + 1
      year = year + 1
    }
    
    yearly_rep$degree_day <- degree_day
    yearly_rep$Tyear_ref <- Tyear_ref
    
    j = Tyear_ref - year0
    temp <- degree_day[j]
    
    yearly_rep$rel_diff <- 1. - yearly_rep$NDVImax / yearly_rep$prediction
    yearly_rep$rel_diffT <- 1. - yearly_rep$degree_day / temp
    
    nNAs <- sum(is.na(yearly_rep$NDVImax))
    if (nNAs < round(0.25 * nyears)) {
      #To be on the safe side, we require for most of the years to give us a non-NA representative or we reject the location.
      
      file <- paste0(folder_out, "NDVIpctCHANGE_site_", site_name, ".png")
      png(file)
      plot(
        yearly_rep$dec_date_ADJ,
        yearly_rep$rel_diff,
        xlab = "Year",
        ylab = "Variation",
        main = "",
        pch = 20,
        type = "l",
        lwd = 2,
        lty = 1,
        col = "blue",
        xaxt = 'n'
      )
      axis(1, at = seq(2000, 2019, by = 1), las = 2)
      dev.off()
      
      file <- paste0(folder_out, "TpctCHANGE_site_", site_name, ".png")
      png(file)
      plot(
        yearly_rep$dec_date_ADJ,
        yearly_rep$rel_diffT,
        xlab = "Year",
        ylab = "Variation",
        main = "",
        pch = 20,
        type = "l",
        lwd = 2,
        lty = 1,
        col = "blue",
        xaxt = 'n'
      )
      axis(1, at = seq(2000, 2019, by = 1), las = 2)
      dev.off()
      
      yearly_rep$date <- as.Date(yearly_rep$date)
      DEF_ALL <- rbind(DEF_ALL, yearly_rep)
    }
    
  } else{
    #IF WATER
    
    j = 1
    yearly_rep <- data.frame(matrix(ncol = length(names(data)), nrow = 0))
    colnames(yearly_rep) <- names(data)
    
    yearly_rep <-
      add_column(yearly_rep, dec_date_ADJ = as.double(), .after = "QA")
    yearly_rep <-
      add_column(yearly_rep, prediction = NA, .after = "dec_date_ADJ")
    yearly_rep <- add_column(yearly_rep, year_ref = NA, .after = "prediction")
    yearly_rep <- add_column(yearly_rep, temp_Kelvin = NA, .after = "year_ref")
    yearly_rep <-
      add_column(yearly_rep, predictionT = NA, .after = "temp_Kelvin")
    yearly_rep <- add_column(yearly_rep, year_refT = NA, .after = "predictionT")
    yearly_rep <-
      add_column(yearly_rep,
                 dec_date_NDVImax = NA,
                 .after = "year_refT")
    yearly_rep <-
      add_column(yearly_rep, NDVImax = NA, .after = "dec_date_NDVImax")
    yearly_rep <- add_column(yearly_rep, degree_day = NA, .after = "NDVImax")
    yearly_rep <- add_column(yearly_rep, Tyear_ref = NA, .after = "degree_day")
    
    year = first_year
    while (year <= last_year) {
      temp_loc <- which(data$year == year)
      yearly_rep[j, seq(1, ncol(data))] <- data[temp_loc[1], ]
      
      yearly_rep$dec_date_ADJ[j] <- year + 0.
      
      j = j + 1
      year = year + 1
    }
    yearly_rep$prediction <- NA
    yearly_rep$year_ref <- NA
    yearly_rep$temp_Kelvin <- NA
    yearly_rep$predictionT <- NA
    yearly_rep$year_refT <- NA
    yearly_rep$dec_date_NDVImax <- NA
    yearly_rep$NDVImax <- NA
    yearly_rep$degree_day <- NA
    yearly_rep$Tyear_ref <- NA
    yearly_rep$rel_diff <- NA
    yearly_rep$rel_diffT <- NA
    
    print("WATER?")
  }
  
  yearly_rep$date <- as.Date(yearly_rep$date)
  DEF_ALL_w_NAs <- rbind(DEF_ALL_w_NAs, yearly_rep)
  
  
}#End locations

#Exporting to a file:

file_out <- paste0(folder_out, "NDVI_and_T-ALL_YEARS.txt")
write.table(
  DEF_ALL,
  file_out,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

#file_out<-paste0(folder_out,"NDVI_and_T-ALL_YEARS_w_NAs.txt")
#write.table(DEF_ALL,file_out, sep="\t",  row.names=FALSE, col.names=TRUE)

#Calculating histograms for NDVI and rel difference, and the associated density functions:

year = first_year
while (year <= last_year) {
  locs <- which(DEF_ALL$year == year)
  temp_temp <- DEF_ALL[locs, ]
  temp <- temp_temp[which(!is.na(temp_temp$NDVImax)), ]
  
  hNDVI <-
    hist(temp$NDVImax,
         plot = FALSE,
         breaks = seq(NDVI_threshold, 1.05, 0.01))
  total <- sum(hNDVI$counts)
  bandwidth <-
    max(hNDVI$breaks) / length(hNDVI$mids)#Even more smoothing, adapted to the histogram plot.
  #dNDVI<-density(temp$NDVImax,bw=bandwidth,na.rm=TRUE)
  ##temp3<-kde(temp$NDVImax)
  ##dNDVI<-data.frame("x"=temp3$eval.points,"y"=temp3$estimate)
  dNDVI <- density(temp$NDVImax, na.rm = TRUE, adjust = 1.25)
  
  NORM <- max(dNDVI$y) / max(hNDVI$counts / total)
  
  file_out = paste0(folder_out, "histogram_NDVI_", year, ".txt")
  OUT = cbind(hNDVI$mids, hNDVI$counts, hNDVI$counts / total, hNDVI$density)
  write.table(
    OUT,
    file_out,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  
  file_out <- paste0(folder_out, "density_NDVI_", year, ".txt")
  OUT = cbind(dNDVI$x, dNDVI$y, -log(dNDVI$y), NORM)
  write.table(
    OUT,
    file_out,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  
  file_out = paste0(folder_out, "histogram_NDVI_", year, ".png")
  df <- data.frame("x" = hNDVI$mids, "y" = hNDVI$counts / total)
  df2 <- data.frame("x" = dNDVI$x, "y" = dNDVI$y / NORM)
  p <-
    ggplot(df, aes(x = x, y = y)) + geom_bar(stat = "identity", fill = "blue") + geom_line(data =
                                                                                             df2, aes(x = x, y = y), colour = "red") + labs(x = "NDVI", y = "Probability") +
    mytheme + scale_x_continuous(limits = c(NDVI_threshold - 0.05, 1.0 + 0.05))
  #scale_x_continuous(limits=c(min(df$x)-0.05,max(df$x)+0.05))
  ggsave(file_out, p, device = "png", dpi = "screen")
  
  relMAX <- max(temp$rel_diff)
  relMIN <- min(temp$rel_diff)
  hrel <-
    hist(
      temp$rel_diff,
      plot = FALSE,
      breaks = seq(relMIN - 0.01, relMAX + 0.01, 0.01)
    )
  total <- sum(hrel$counts)
  bandwidth <-
    max(hrel$breaks) / length(hrel$mids)#Even more smoothing, adapted to the histogram plot.
  #drel<-density(temp$rel_diff,bw=bandwidth,na.rm=TRUE)
  ##temp3<-kde(temp$rel_diff)
  ##drel<-data.frame("x"=temp3$eval.points,"y"=temp3$estimate)
  drel <- density(temp$rel_diff, na.rm = TRUE, adjust = 1.25)
  
  NORM <- max(drel$y) / max(hrel$counts / total)
  
  file_out = paste0(folder_out, "histogram_rel_diff_", year, ".txt")
  OUT = cbind(hrel$mids, hrel$counts, hrel$counts / total, hrel$density)
  write.table(
    OUT,
    file_out,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  
  file_out <- paste0(folder_out, "density_rel_diff_", year, ".txt")
  OUT = cbind(drel$x, drel$y, -log(drel$y), NORM)
  write.table(
    OUT,
    file_out,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  
  file_out = paste0(folder_out, "histogram_rel_diff_", year, ".png")
  df <- data.frame("x" = hrel$mids, "y" = hrel$counts / total)
  df2 <- data.frame("x" = drel$x, "y" = drel$y / NORM)
  p <-
    ggplot(df, aes(x = x, y = y)) + geom_bar(stat = "identity", fill = "blue") + geom_line(data =
                                                                                             df2, aes(x = x, y = y), colour = "red") + labs(x = "Relative difference", y =
                                                                                                                                              "Probability") + mytheme + scale_x_continuous(limits = c(-0.5 - 0.05, 0.5 +
                                                                                                                                                                                                         0.05))
  #difference",y="Probability")+mytheme+ scale_x_continuous(limits=c(min(df$x)-0.05,max(df$x)+0.05))
  ggsave(file_out, p, device = "png", dpi = "screen")
  
  temp2 <- unlist(lapply(temp$rel_diff, step_f))
  hrelPOS <-
    hist(
      temp2,
      plot = FALSE,
      breaks = seq(0.0, max(temp2) + 0.05, 0.01),
      na.rm = TRUE
    )
  total <- sum(hrelPOS$counts)
  bandwidth <-
    max(hrelPOS$breaks) / length(hrelPOS$mids)#Even more smoothing, adapted to the histogram plot.
  #drelPOS<-density(temp2,bw=bandwidth,na.rm=TRUE)
  ##temp3<-kde(temp2)
  ##drelPOS<-data.frame("x"=temp3$eval.points,"y"=temp3$estimate)
  ##temp3<-logspline(temp2,lbound = 0)
  drelPOS <- density(temp2,
                     adjust = 1.25,
                     from = 0,
                     na.rm = TRUE)
  ##plot(temp3,  col = "red",  lwd = 3, add = TRUE)
  
  NORM <- max(drelPOS$y) / max(hrelPOS$counts / total)
  
  file_out = paste0(folder_out, "histogram_rel_diff_LUMPED_", year, ".txt")
  OUT = cbind(hrelPOS$mids,
              hrelPOS$counts,
              hrelPOS$counts / total,
              hrelPOS$density)
  write.table(
    OUT,
    file_out,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  
  file_out <- paste0(folder_out, "density_rel_diff_LUMPED_", year, ".txt")
  OUT = cbind(drelPOS$x, drelPOS$y, -log(drelPOS$y), NORM)
  write.table(
    OUT,
    file_out,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  
  file_out = paste0(folder_out, "histogram_rel_diff_LUMPED_", year, ".png")
  df <- data.frame("x" = hrelPOS$mids, "y" = hrelPOS$counts / total)
  df2 <- data.frame("x" = drelPOS$x, "y" = drelPOS$y / NORM)
  p <-
    ggplot(df, aes(x = x, y = y)) + geom_bar(stat = "identity", fill = "blue") + geom_line(data =
                                                                                             df2, aes(x = x, y = y), colour = "red") + labs(x = "Relative difference", y =
                                                                                                                                              "Probability") + mytheme + scale_x_continuous(limits = c(0.0 - 0.05, 0.5 +
                                                                                                                                                                                                         0.05))
  #difference",y="Probability")+mytheme+ scale_x_continuous(limits=c(min(df$x)-0.05,max(df$x)+0.05))
  ggsave(file_out, p, device = "png", dpi = "screen")
  
  
  year = year + 1
  
}#End of years


#Histogram with all years of reference (the hope is that it'll be very narrowly distributed):
#Because it's the same number for all years of a location, it's faster to pick a year than to sweep all locations:

year = first_year + floor(nyears / 2)

locs <- which(DEF_ALL$year == year)
temp_temp <- DEF_ALL[locs, ]
temp <- temp_temp[which(!is.na(temp_temp$year_ref)), ]

hyear <-
  hist(temp$year_ref,
       plot = FALSE,
       breaks = seq(first_year - 1, last_year + 1, 1))
total <- sum(hyear$counts)
bandwidth <-
  max(hyear$breaks) / length(hyear$mids)#Even more smoothing, adapted to the histogram plot.
#dyear<-density(temp$year_ref,bw=bandwidth,na.rm=TRUE)
##temp3<-kde(temp$year_ref)
##dyear<-data.frame("x"=temp3$eval.points,"y"=temp3$estimate)
dyear <- density(temp$year_ref, na.rm = TRUE, adjust = 1.25)

NORM <- max(dyear$y) / max(hyear$counts / total)

file_out = paste0(folder_out, "histogram_ref_year.txt")
OUT = cbind(hyear$mids, hyear$counts, hyear$counts / total, hyear$density)
write.table(
  OUT,
  file_out,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

file_out <- paste0(folder_out, "density_ref_year.txt")
OUT = cbind(dyear$x, dyear$y, -log(dyear$y), NORM)
write.table(
  OUT,
  file_out,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

file_out = paste0(folder_out, "histogram_ref_year.png")
df <- data.frame("x" = hyear$mids, "y" = hyear$counts / total)
df2 <- data.frame("x" = dyear$x, "y" = dyear$y / NORM)
p <-
  ggplot(df, aes(x = x, y = y)) + geom_bar(stat = "identity", fill = "blue") + geom_line(data =
                                                                                           df2, aes(x = x, y = y), colour = "red") + labs(x = "Reference year", y =
                                                                                                                                            "Probability") + mytheme + scale_x_continuous(limits = range(hyear$mids))
ggsave(file_out, p, device = "png", dpi = "screen")



#PLOT THE MAPS FOR EACH YEAR: https://slcladal.github.io/maps.html

#Adding the RGB color column. First, the following function rbPal(x) will generate a gradient of x=ncolors colors between the two colors you provide:
rbPal <- colorRampPalette(c('darkolivegreen', 'darkred'))
ncolors <-
  11#So that we can do the "fixed/normalized" color palette more easily, as a number between 0 and 10.

#Then, assign each rel_diff a value from 1 (minimum rel_diff) to ncolors (max rel_diff) and assign one of the ncolors colors. We choose green if no defoliation (rel_diff<=0) and red if there is defoliation:

#If there are unreasonably high values of rel_diff (rel_diff>1) due to bad NDVI measure/month selection/model prediction, they will screw the color code/palette. In those cases we make sure rel_diff=1 is the maximum:

DEF_ALL_w_NAs$col <-
  rbPal(ncolors)[as.numeric(cut(unlist(lapply(
    DEF_ALL_w_NAs$rel_diff, step_f2
  )), breaks = ncolors))]

#Then, we assing blue to NaNs, and black to unreasonable values (e.g. rel_diff > 1.):
DEF_ALL_w_NAs$col[which(is.na(DEF_ALL_w_NAs$rel_diff))] <-
  col2hex("blue")
DEF_ALL_w_NAs$col[which(DEF_ALL_w_NAs$rel_diff > 1)] <-
  col2hex("black")

#Now we produce a version of the color palette that is fixed between 0 and MAXrel:
####MAXrel<-max(unlist(lapply(DEF_ALL_w_NAs$rel_diff,step_f)))
####DEF_ALL_w_NAs$col_norm <- rbPal(ncolors)[as.numeric(cut(unlist(lapply(DEF_ALL_w_NAs$rel_diff,step_f))/MAXrel,breaks = ncolors))]#Nope, that would be like you had it before, with the max being the last color.

MAXrel <- 0.5
DEF_ALL_w_NAs$col_norm <-
  rbPal(ncolors)[1 + round(10 * unlist(lapply(DEF_ALL_w_NAs$rel_diff / MAXrel, step_f2)))]
DEF_ALL_w_NAs$col_norm[which(is.na(DEF_ALL_w_NAs$rel_diff))] <-
  col2hex("blue")
DEF_ALL_w_NAs$col_norm[which(DEF_ALL_w_NAs$rel_diff > 1)] <-
  col2hex("black")

year = year0
while (year <= last_year) {
  locs <- which(DEF_ALL_w_NAs$year == year)
  
  bbox <-
    make_bbox(
      lon = c(long_min - 0.1, long_max + 0.1),
      lat = c(lat_min - 0.05, lat_max + 0.05),
      f = .1
    )
  region <-
    get_map(
      location = bbox,
      zoom = 10,
      maptype = "terrain",
      style = 'feature:administrative.country|element:labels|visibility:off',
      legend = "topleft"
    )
  #region<-get_map(location=bbox, zoom=5, maptype="terrain", color="bw")
  #region<-get_map(location=bbox, zoom=10, maptype="terrain", color="bw")
  
  #p<-ggmap(region)+geom_point(data = DEF_ALL_w_NAs[locs,], mapping = aes(x = DEF_ALL_w_NAs$long_dd[locs], y = DEF_ALL_w_NAs$lat_dd[locs]), color = DEF_ALL_w_NAs$col[locs],shape=15, size=3)+labs(x="Longitude",y="Latitude")
  #FOR A PLOT IN WHICH YOU DON'T HAVE GAPS BETWEEN LOCATIONS (BUT IT'S MUCH SLOWER TO PLOT!!):
  #p<-ggmap(region)+geom_tile(data = DEF_ALL_w_NAs[locs,], mapping = aes(x = DEF_ALL_w_NAs$long_dd[locs], y = DEF_ALL_w_NAs$lat_dd[locs]), fill=DEF_ALL_w_NAs$col[locs])+labs(x="Longitude",y="Latitude")
  
  #EASIER VERSION, WITHOUT THE NEED OF THE COLOR COLUMN (if you invoke the library "scales" and add ",oob=squish)" the squish part makes the points that are out of bounds look like the limits):
  
  p <-
    ggmap(region) + geom_point(
      data = DEF_ALL_w_NAs[locs, ],
      mapping = aes(
        x = DEF_ALL_w_NAs$long_dd[locs],
        y = DEF_ALL_w_NAs$lat_dd[locs],
        color = unlist(lapply(DEF_ALL_w_NAs$rel_diff[locs], step_f3)) + 0.001
      ),
      shape = 15,
      size = 3
    ) + labs(x = "Longitude", y = "Latitude", color = "rel_diff") + scale_color_gradient(
      low = "green",
      high = "darkred",
      na.value = "black"
    )
  
  #SAVE TO FILE: https://ggplot2.tidyverse.org/reference/ggsave.html
  file <- paste0(folder_out, "map_NDVIpctCHANGE_", year, ".png")
  ggsave(file, p, device = "png", dpi = "screen")
  
  #p<-ggmap(region)+geom_point(data = DEF_ALL_w_NAs[locs,], mapping = aes(x = DEF_ALL_w_NAs$long_dd[locs], y = DEF_ALL_w_NAs$lat_dd[locs]), color = DEF_ALL_w_NAs$col_norm[locs],shape=15, size=3)+labs(x="Longitude",y="Latitude")
  #FOR A PLOT IN WHICH YOU DON'T HAVE GAPS BETWEEN LOCATIONS (BUT IT'S MUCH SLOWER TO PLOT!!):
  #p<-ggmap(region)+geom_tile(data = DEF_ALL_w_NAs[locs,], mapping = aes(x = DEF_ALL_w_NAs$long_dd[locs], y = DEF_ALL_w_NAs$lat_dd[locs]), fill=DEF_ALL_w_NAs$col_norm[locs])+labs(x="Longitude",y="Latitude")
  
  MAXrel <- 0.5
  p <-
    ggmap(region) + geom_point(
      data = DEF_ALL_w_NAs[locs, ],
      mapping = aes(
        x = DEF_ALL_w_NAs$long_dd[locs],
        y = DEF_ALL_w_NAs$lat_dd[locs],
        color = unlist(lapply(
          DEF_ALL_w_NAs$rel_diff[locs] / MAXrel, step_f3
        )) + 0.001
      ),
      shape = 15,
      size = 3
    ) + labs(x = "Longitude", y = "Latitude", color = "rel_diff") + scale_color_gradient(
      low = "green",
      high = "darkred",
      na.value = "black",
      limits = c(0, 1)
    )
  
  #SAVE TO FILE: https://ggplot2.tidyverse.org/reference/ggsave.html
  file <- paste0(folder_out, "map_NDVIpctCHANGE_FIX_", year, ".png")
  ggsave(file, p, device = "png", dpi = "screen")
  
  #DEGREE-DAY MAPS!!!
  
  p <-
    ggmap(region) + geom_point(
      data = DEF_ALL_w_NAs[locs, ],
      mapping = aes(
        x = DEF_ALL_w_NAs$long_dd[locs],
        y = DEF_ALL_w_NAs$lat_dd[locs],
        color = DEF_ALL_w_NAs$rel_diffT[locs]
      ),
      shape = 15,
      size = 3
    ) + labs(x = "Longitude", y = "Latitude", color = "rel_diff") + scale_color_gradient(low =
                                                                                           "red",
                                                                                         high = "blue",
                                                                                         na.value = "black")
  
  #SAVE TO FILE: https://ggplot2.tidyverse.org/reference/ggsave.html
  file <- paste0(folder_out, "map_DEGREEDAYpctCHANGE_", year, ".png")
  ggsave(file, p, device = "png", dpi = "screen")
  
  p <-
    ggmap(region) + geom_point(
      data = DEF_ALL_w_NAs[locs, ],
      mapping = aes(
        x = DEF_ALL_w_NAs$long_dd[locs],
        y = DEF_ALL_w_NAs$lat_dd[locs],
        color = DEF_ALL_w_NAs$rel_diffT[locs]
      ),
      shape = 15,
      size = 3
    ) + labs(x = "Longitude", y = "Latitude", color = "rel_diff") + scale_color_gradient(
      low = "red",
      high = "blue",
      na.value = "black",
      limits = c(-2, 2)
    )
  
  #SAVE TO FILE: https://ggplot2.tidyverse.org/reference/ggsave.html
  file <- paste0(folder_out, "map_DEGREEDAYpctCHANGE_FIX_", year, ".png")
  ggsave(file, p, device = "png", dpi = "screen")
  
  
  year = year + 1
}



###############CLASSIFICATION OF DEFOLIATION LEVEL (ORDER PARAMETER, p) BY VALUE OF THE ACCUMULATED DEGREE-DAY (CONTROL PARAMETER, r)################################
###################################REMEMBER: p > 0 MEANS DEFOLIATION, AND r < 0 MEANS HOTTER THAN THE REFERENCE######################################################


rmin <-
  floor(min(na.omit(DEF_ALL$rel_diffT)))    #All processed dataframes have the same final set of climate variables, so we can use this binning for any of them.
rmax <- ceiling(max(na.omit(DEF_ALL$rel_diffT)))
drs <- c(0.1, 0.2, 0.3)

for (dr in drs) {
  addfolder = paste0("dr_", dr, "_NOT_MERGING_nmin200/")
  dir.create(file.path(mainDir, folder_out, addfolder))
  
  nbins <- ceiling((rmax - rmin) / dr) + 2
  breaksr <- seq(from = rmin - dr,
                 by = dr,
                 length.out = nbins)
  temp <- tail(breaksr, n = 1)
  if (temp < rmax + dr) {
    breaksr <-
      unique(append(breaksr, seq(temp, rmax + dr, dr)))#So the bins will cover from breaksr to breaksr+dr (the latter excluded)
  }
  middlers <- 0
  for (i in 1:(length(breaksr) - 1)) {
    middlers[i] <- (breaksr[i] + breaksr[i + 1]) / 2
  }
  middlers[i + 1] <- breaksr[i + 1] + dr / 2
  
  samer <- data.frame(matrix(ncol = length(breaksr), nrow = 0))
  
  leftr <-
    round(middlers - dr / 2, digits = 10) #To avoid quasi-zeros and similar
  rightr <-
    round(middlers + dr / 2, digits = 10) #To avoid quasi-zeros and similar
  
  colnames(samer) <-
    paste0("r_", trimws(format(leftr, nsmall = 2)), "-", trimws(format(rightr, nsmall =
                                                                         2)))
  
  for (j in 1:nrow(DEF_ALL))
  {
    ir <-
      .bincode(DEF_ALL$rel_diffT[j],
               breaksr,
               right = TRUE,
               include.lowest = FALSE)
    samer[sum(!is.na(samer[, ir])) + 1, ir] <- DEF_ALL$rel_diff[j]
    
    print(paste0("dr=", dr, " -> ", j, " out of ", nrow(DEF_ALL)))
    
  }
  
  nmin = 200
  
  #keep a list of the bins that will ultimately be removed or merged:
  list_removed <- colSums(is.na(samer)) >= nrow(samer) - nmin
  
  #samer<- samer[,colSums(is.na(samer))<nrow(samer)-nmin]#We keep only columns with at least nmin+1 points
  samer <-
    samer[, colSums(is.na(samer)) < nrow(samer)]#We keep only columns with at least 1 point
  
  
  
  #MERGING: Expecting that only extremes r-bins will show low counts, we merge those with the first bin next to them with sufficient counts:
  #bin_low<-min(which(colSums(is.na(samer))<nrow(samer)-nmin))
  #bin_high<-max(which(colSums(is.na(samer))<nrow(samer)-nmin))
  #bin_max<-ncol(samer)
  #
  #for (j in 1:(bin_low-1))
  #{
  #for (i in 1:sum(!is.na(samer[,j]))){
  #
  #samer[sum(!is.na(samer[,bin_low]))+1,bin_low]<-samer[i,j]
  ##print(paste(i,j))
  #}
  #}
  #
  #for (j in (bin_high+1):bin_max)
  #{
  #for (i in 1:sum(!is.na(samer[,j]))){
  #
  #samer[sum(!is.na(samer[,bin_high]))+1,bin_high]<-samer[i,j]
  ##print(paste(i,j))
  #}
  #}
  
  
  
  #And remove the bins that are now "empty":
  samer <-
    samer[, colSums(is.na(samer)) < nrow(samer) - nmin]#We keep only columns with at least nmin+1 points
  
  #Calculate and export the potentials:
  
  
  ncolors <- length(colnames(samer))
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  colors <- sample(color, ncolors)
  #rbPal <- colorRampPalette(c('blue','red'))
  #colors <- rbPal(ncolors)
  
  
  file <- paste0(folder_out, addfolder, "ALL_defoliation_densities.png")
  png(file)
  plot.new()
  plot.window(xlim = c(-0.6, 0.6), ylim = c(0, 0.14))
  axis(1)
  axis(2)
  box()
  
  for (i in colnames(samer)) {
    temp <- samer[, i]
    column <- samer[!is.na(temp), i]
    relMAX <- max(column)
    relMIN <- min(column)
    H <-
      hist(column,
           plot = FALSE,
           breaks = seq(relMIN - 0.01, relMAX + 0.01, 0.01))
    D <- density(column, na.rm = T)
    #D<-density(column,na.rm=T,adjust=1.25)
    #temp3<-kde(column)
    #D<-data.frame("x"=temp3$eval.points,"y"=temp3$estimate)
    total <- sum(H$counts)
    
    NORM <- max(D$y) / max(H$counts / total)
    
    file_out <- paste0(folder_out, addfolder, "histogram_", i, ".txt")
    OUT = cbind(H$mids, H$counts, H$counts / total, H$density, total)
    write.table(
      OUT,
      file_out,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
    file_out <- paste0(folder_out, addfolder, "density_", i, ".txt")
    OUT = cbind(D$x, D$y, -log(D$y), NORM)
    write.table(
      OUT,
      file_out,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
    
    #Plot histograms and density curves:
    
    file_out = paste0(folder_out, addfolder, "histogram_", i, ".png")
    df <- data.frame("x" = H$mids, "y" = H$counts / total)
    df2 <- data.frame("x" = D$x, "y" = D$y / NORM)
    p <-
      ggplot(df, aes(x = x, y = y)) + geom_bar(stat = "identity", fill = "blue") + geom_line(data =
                                                                                               df2, aes(x = x, y = y), colour = "red") + labs(x = "defoliation (%)", y =
                                                                                                                                                "Probability") + mytheme + scale_x_continuous(limits = c(-0.6, 0.6)) + scale_y_continuous(limits =
                                                                                                                                                                                                                                            c(0., 0.14))
    ggsave(file_out, p, device = "png", dpi = "screen")
    
    #Plot all density curves together:
    j = which(colnames(samer) == i)
    #lines(D$x,D$y/NORM,xlab="defoliation (%)", ylab="Probability",col=sample(color,1))
    #lines(D$x,D$y/NORM,xlab="defoliation (%)", ylab="Probability",col=sample(color,1))
    lines(D$x, D$y / NORM, col = colors[j], lwd = 2)
  }
  #legend("topright", legend = colnames(samer),col=colors,lwd=2,ncol=2)
  legend(
    "topright",
    legend = colnames(samer),
    col = colors,
    lwd = 1,
    cex = 0.6
  )
  text(
    x = 0,
    y = -0.03,
    labels = "defoliation (%)",
    xpd = NA
  )
  text(
    x = -0.6 - 0.2,
    y = 0.07,
    labels = "Probability",
    srt = 90,
    xpd = NA
  )
  dev.off()
  
  #Sketch the phase diagram:
  file <- paste0(folder_out, addfolder, "phase_diagram_dr_", dr, ".png")
  png(file)
  nremoved <- length(middlers[which(list_removed == FALSE)])
  rmin_eff <- middlers[which(list_removed == FALSE)[[1]]] - 0.01
  rmax_eff <- middlers[which(list_removed == FALSE)[[nremoved]]] + 0.01
  
  j = 1
  k = 1
  phase_x <- 0
  phase_y <- 0
  for (i in colnames(samer)) {
    temp <- samer[, i]
    column <- samer[!is.na(temp), i]
    ycoord <- Modes(column)$modes
    xcoord <- middlers[which(list_removed == FALSE)[[j]]]
    phase_x[k] <- xcoord
    phase_y[k] <- ycoord[1]
    k = k + 1
    if (length(ycoord) > 1) {
      #If there's bimodality!!!! (check the documentation for the Modes function, because you can play with the smoothing and criteria to reject modes)
      phase_x[k] <- xcoord
      phase_y[k] <- ycoord[2]
      k = k + 1
    }
    print(paste(xcoord, ycoord))
    j = j + 1
  }
  plot(
    phase_x,
    phase_y,
    xlim = c(rmin_eff, rmax_eff),
    ylim = c(0, 0.15),
    xlab = "climatic factor, r",
    ylab = "Mode of defoliation (%)"
  )
  dev.off()
  
  
  smooth <- predict(loess(phase_y ~ phase_x))
  smooth2 <-
    predict(loess.as(
      phase_x,
      phase_y,
      degree = 1,
      criterion = c("aicc", "gcv")[2],
      user.span = NULL,
      plot = F
    ))
  file_out <- paste0(folder_out, addfolder, "phase_diagram_dr_", dr, ".txt")
  OUT = cbind(phase_x, phase_y, smooth, smooth2)
  write.table(
    OUT,
    file_out,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  
}
