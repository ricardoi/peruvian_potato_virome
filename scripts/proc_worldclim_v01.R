#'@ authors: Alcala-Briseno, Andersen-Onofre and Garrett 
#'@ title: Processing world clim data
#'@ WorldClim version :2.1 climate data for 1970-2000
#'@ released: Jan 2020
#'@ version: Rstudio v.4.0.3

# set working directory
setwd("~/Dropbox (UFL)/Alcala_Briseno-Garrett/++Papa_virome/papa/")

# load libraries
library(raster)
library(tidyverse)
library(usdm)

# monthly climate data
# 19 GeoTiff (.tif)
# time: 30 seconds (~1 km2) ]
# wordclim version: 2.1

# variables:
{# minimum temperature (°C)	      | tmin |
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
}

# load incidence matrix 
# potato virome presence absence data
ppv <- as_tibble(read.csv("papa_Data_Virome_v2.csv")) %>%
        dplyr::select(Latitude, Longitude, 4)

# isolate latitude and longitude for each sample locations
head(ppv)
# lat lot virus.sp

# note: take into account the virome nature of the data

# reading worldclim data 
vars <- c("tmin", "tmax", "tavg", "prec", "srad", "wind", "vapr")
bioclim <- paste0("BIO", 1:19)

# creating lists
# for bioclimatic variables and elevation
tmin = tmax = tavg = prec = srad = wind = vapr = bio = elev = list()

# tmin
dirs <- list.dirs("6-spatial_analysis/WorldClim_v2.1/wc2.1_30s_tmin/")
tmins = list.files(dirs, pattern="*.tif", full.names = T, all.files = T)
for (t in seq_along(tmins)){
  tmin[[t]] <- raster(tmins[[t]])
}
tmins <- raster::stack(tmin[1], tmin[2], tmin[3], tmin[4], tmin[5], tmin[6], 
                       tmin[7], tmin[8], tmin[9], tmin[10], tmin[11], tmin[12])

# tmax
dirs <- list.dirs("6-spatial_analysis/WorldClim_v2.1/wc2.1_30s_tmax/")
tmaxs = list.files(dirs, pattern="*.tif", full.names = T, all.files = T)
for (t in seq_along(tmaxs)){
  tmax[[t]] <- raster(tmaxs[[t]])
}
tmaxs <- raster::stack(tmax[1], tmax[2], tmax[3], tmax[4], tmax[5], tmax[6], 
                       tmax[7], tmax[8], tmax[9], tmax[10], tmax[11], tmax[12])
# tavg
dirs <- list.dirs("6-spatial_analysis/WorldClim_v2.1/wc2.1_30s_tavg/")
tavgs = list.files(dirs, pattern="*.tif", full.names = T, all.files = T)
for (t in seq_along(tavgs)){
  tavg[[t]] <- raster(tavgs[[t]])
}
tavgs <- raster::stack(tavg[1], tavg[2], tavg[3], tavg[4], tavg[5], tavg[6], 
                       tavg[7], tavg[8], tavg[9], tavg[10], tavg[11], tavg[12])
# prec
dirs <- list.dirs("6-spatial_analysis/WorldClim_v2.1/wc2.1_30s_prec/")
precs = list.files(dirs, pattern="*.tif", full.names = T, all.files = T)
for (t in seq_along(precs)){
  prec[[t]] <- raster(precs[[t]])
}
precs <- raster::stack(prec[1], prec[2], prec[3], prec[4], prec[5], prec[6], 
                       prec[7], prec[8], prec[9], prec[10], prec[11], prec[12])
# srad
dirs <- list.dirs("6-spatial_analysis/WorldClim_v2.1/wc2.1_30s_srad/")
srads = list.files(dirs, pattern="*.tif", full.names = T, all.files = T)
for (t in seq_along(srads)){
  srad[[t]] <- raster(srads[[t]])
}
srads <- raster::stack(srad[1], srad[2], srad[3], srad[4], srad[5], srad[6], 
                       srad[7], srad[8], srad[9], srad[10], srad[11], srad[12])
# wind
dirs <- list.dirs("6-spatial_analysis/WorldClim_v2.1/wc2.1_30s_wind/")
winds = list.files(dirs, pattern="*.tif", full.names = T, all.files = T)
for (t in seq_along(winds)){
  wind[[t]] <- raster(winds[[t]])
}
winds <- raster::stack(wind[1], wind[2], wind[3], wind[4], wind[5], wind[6], 
                       wind[7], wind[8], wind[9], wind[10], wind[11], wind[12])
# vapr
dirs <- list.dirs("6-spatial_analysis/WorldClim_v2.1/wc2.1_30s_vapr/")
vaprs = list.files(dirs, pattern="*.tif", full.names = T, all.files = T)
for (t in seq_along(vaprs)){
  vapr[[t]] <- raster(vaprs[[t]])
}
vaprs <- raster::stack(vapr[1], vapr[2], vapr[3], vapr[4], vapr[5], vapr[6], 
                       vapr[7], vapr[8], vapr[9], vapr[10], vapr[11], vapr[12])
#---- BIOCLIM
# bio
dirs <- list.dirs("6-spatial_analysis/WorldClim_v2.1/wc2.1_30s_bio/")
bios = list.files(dirs, pattern="*.tif", full.names = T, all.files = T)
for (t in seq_along(bios)){
  bio[[t]] <- raster(bios[[t]])
}
bios <- raster::stack(bio[1], bio[2], bio[3], bio[4], bio[5], bio[6], 
                      bio[7], bio[8], bio[9], bio[10], bio[11], bio[12],
                      bio[13], bio[14], bio[15], bio[16], bio[17], bio[18],
                      bio[19])
# elev
dirs <- list.dirs("6-spatial_analysis/WorldClim_v2.1/wc2.1_30s_elev/")
elevs = list.files(dirs, pattern="*.tif", full.names = T, all.files = T)
for (t in seq_along(elevs)){
  elev[[t]] <- raster(elevs[[t]])
}
elevs <- raster::stack(elev[1])

### then stack all of the stacks 

all.vars <- raster::stack(tmin[[1]], tmin[[2]], tmin[[3]], tmin[[4]], tmin[[5]], tmin[[6]],
                          tmin[[7]], tmin[[8]], tmin[[9]], tmin[[10]], tmin[[11]], tmin[[12]],
                          tmax[[1]], tmax[[2]], tmax[[3]], tmax[[4]], tmax[[5]], tmax[[6]],
                          tmax[[7]], tmax[[8]], tmax[[9]], tmax[[10]], tmax[[11]], tmax[[12]],
                          tavg[[1]], tavg[[2]], tavg[[3]], tavg[[4]], tavg[[5]], tavg[[6]],
                          tavg[[7]], tavg[[8]], tavg[[9]], tavg[[10]], tavg[[11]], tavg[[12]],
                          prec[[1]], prec[[2]], prec[[3]], prec[[4]], prec[[5]], prec[[6]],
                          prec[[7]], prec[[8]], prec[[9]], prec[[10]], prec[[11]], prec[[12]],
                          srad[[1]], srad[[2]], srad[[3]], srad[[4]], srad[[5]], srad[[6]],
                          srad[[7]], srad[[8]], srad[[9]], srad[[10]], srad[[11]], srad[[12]],
                          wind[[1]], wind[[2]], wind[[3]], wind[[4]], wind[[5]], wind[[6]],
                          wind[[7]], wind[[8]], wind[[9]], wind[[10]], wind[[11]], wind[[12]],
                          vapr[[1]], vapr[[2]], vapr[[3]], vapr[[4]], vapr[[5]], vapr[[6]],
                          vapr[[7]], vapr[[8]], vapr[[9]], vapr[[10]], vapr[[11]], vapr[[12]],
                          bio[[1]], bio[[2]], bio[[3]], bio[[4]], bio[[5]], bio[[6]],
                          bio[[7]], bio[[8]], bio[[9]], bio[[10]], bio[[11]], bio[[12]],
                          elev[[1]])




# select data based on disease incidence 
extractedVIF <- raster::extract(all.vars, ppv[,1:2])
# 
vifcor(extractedVIF,th=0.95)

# remove hightly coorelated variables
VIFvars <- vifstep(extractedVIF, th = 60)

### now create a new stack with just the reduced set of predictors 

#Could be something like: 

condensedStacka <- raster::stack(prec[[1]], # wc2.1_30s_prec_01 
                                 prec[[8]], #wc2.1_30s_prec_08
                                 prec[[11]], #wc2.1_30s_prec_11
                                 prec[[12]], #wc2.1_30s_prec_12
                                 srad[[1]], #wc2.1_30s_srad_01
                                 srad[[2]], #wc2.1_30s_srad_02
                                 srad[[3]], #wc2.1_30s_srad_03
                                 srad[[4]], #wc2.1_30s_srad_04
                                 srad[[6]], #wc2.1_30s_srad_06
                                 srad[[7]], #wc2.1_30s_srad_07
                                 srad[[9]], #wc2.1_30s_srad_09
                                 srad[[10]], #wc2.1_30s_srad_10
                                 srad[[12]], #wc2.1_30s_srad_12
                                 vapr[[3]], #wc2.1_30s_vapr_03
                                 vapr[[5]], #wc2.1_30s_vapr_05
                                 vapr[[7]], #wc2.1_30s_vapr_07
                                 vapr[[8]], #wc2.1_30s_vapr_08
                                 vapr[[9]], #wc2.1_30s_vapr_09
                                 vapr[[10]], #wc2.1_30s_vapr_10
                                 bio[[14]], #wc2.1_30s_bio_14
                                 bio[[2]]) #wc2.1_30s_bio_2)
