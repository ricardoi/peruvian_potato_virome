#'@ authors: Alcala-Briseno, Andersen-Onofre and Garrett 
#'@ title: Processing world clim data
#'@ WorldClim version :2.1 climate data for 1970-2000
#'@ released: Jan 2020
#'@ version: Rstudio v.4.0.3

# set working directory
setwd("~/Dropbox (UFL)/Alcala_Briseno-Garrett/++Papa_virome/papa/")

# load libraries
library(raster)


# monthly climate data
# 19 GeoTiff (.tif)
# time: 30 seconds (~1 km2) ]
# wordclim version: 2.1

# variables:
# minimum temperature (°C)	      | tmin |
# maximum temperature (°C)	      | tmax |
# average temperature (°C)	      | tavg |
# precipitation (mm)              | prec |
# solar radiation (kJ m-2 day-1)	| srad |
# wind speed (m s-1)              |	wind |
# water vapor pressure (kPa)	    | vapr |

# Bioclimatic variables
# standard (19) WorldClim Bioclimatic variables for WorldClim version 2.
# average for the years 1970-2000

# They are coded as follows:
# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (×100)
# BIO4 = Temperature Seasonality (standard deviation ×100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter

# load incidence matrix 
# potato virome presence absence data
ppv <- read.csv("peruvian_potato_virome.csv")

# isolate latitude and longitude for each sample locations
head(ppv)
# lat lot virus.sp

# note: take into account the virome nature of the data

# reading worldclim data 
vars <- c("tmin", "tmax", "tavg", "prec", "srad", "wind", "vapr")
bioclim <- paste0("BIO", 1:19)

tmin = tmax = tavg = prec = srad = wind = vapr = list()

# tmin


wc2.0_30s_tmin_01<- raster("Data5/wc2.0_30s_srad/wc2.0_30s_tmin_01.tif")
wc2.0_30s_srad_02<- raster("Data5/wc2.0_30s_srad/wc2.0_30s_srad_02.tif")
wc2.0_30s_srad_03<- raster("Data5/wc2.0_30s_srad/wc2.0_30s_srad_03.tif")
wc2.0_30s_srad_04<- raster("Data5/wc2.0_30s_srad/wc2.0_30s_srad_04.tif")
wc2.0_30s_srad_05<- raster("Data5/wc2.0_30s_srad/wc2.0_30s_srad_05.tif")
wc2.0_30s_srad_06<- raster("Data5/wc2.0_30s_srad/wc2.0_30s_srad_06.tif")
wc2.0_30s_srad_07<- raster("Data5/wc2.0_30s_srad/wc2.0_30s_srad_07.tif")
wc2.0_30s_srad_08<- raster("Data5/wc2.0_30s_srad/wc2.0_30s_srad_08.tif")
wc2.0_30s_srad_09<- raster("Data5/wc2.0_30s_srad/wc2.0_30s_srad_09.tif")
wc2.0_30s_srad_10<- raster("Data5/wc2.0_30s_srad/wc2.0_30s_srad_10.tif")
wc2.0_30s_srad_11<- raster("Data5/wc2.0_30s_srad/wc2.0_30s_srad_11.tif")
wc2.0_30s_srad_12<- raster("Data5/wc2.0_30s_srad/wc2.0_30s_srad_12.tif")


# stack rasters
sradStack <- raster::stack(wc2.0_30s_srad_01, wc2.0_30s_srad_02, wc2.0_30s_srad_03, wc2.0_30s_srad_04, wc2.0_30s_srad_05, wc2.0_30s_srad_06, wc2.0_30s_srad_07, wc2.0_30s_srad_08, wc2.0_30s_srad_09, wc2.0_30s_srad_10, wc2.0_30s_srad_11, wc2.0_30s_srad_12)


