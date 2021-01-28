#'@ authors: Alcala-Briseno, Andersen-Onofre and Garrett 
#'@ title: Processing world clim data
#'@ WorldClim version :2.1 climate data for 1970-2000
#'@ released: Jan 2020
#'@ version: Rstudio v.4.0.3

# set working directory
setwd("~/Dropbox (UFL)/Alcala_Briseno-Garrett/++Papa_virome/papa/")

# monthly climate data
# 12 GeoTiff (.tif)
# 30 seconds (~1 km2) to 10 minutes (~340 km2).

# variables:
# minimum temperature (°C)	      | tmin |
# maximum temperature (°C)	      | tmax |
# average temperature (°C)	      | tavg |
# precipitation (mm)              | prec |
# solar radiation (kJ m-2 day-1)	| srad |
# wind speed (m s-1)              |	wind |
# water vapor pressure (kPa)	    | vapr |
# time: 10m	tmin 5m	tmin 2.5m	tmin 30s
# wordclim version: 2.1

vars <- c("tmin", "tmax", "tavg", "prec", "srad", "wind", "vapr")
time <- "30s"
wc.ver <- "2.1"

# get urls
urls = dests = list()
for (i in seq_along(vars)){
urls[i] <-  paste0("http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc", wc.ver, "_",
       time, "_", vars[i], ".zip")
dests[i] <- paste0("6-spatial_analysis/wc", wc.ver, "_",
                   time, "_", vars[i], ".zip") 
}

for (j in seq_along(vars)){
download.file(urls[[j]], dests[[j]])
}

for (k in seq_along(vars)){
  unzip(dests[[k]])
}

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

url <- "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_bio.zip"
dest <- "6-spatial_analysis/wc2.1_30s_bio.zip"
download.file(url, dest)
unzip("6-spatial_analysis/wc2.1_30s_bio.zip")

# elevation

# derived from the SRTM elevation data.
url <- "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_elev.zip"
dest <- "6-spatial_analysis/wc2.1_30s_elev.zip"
unzip("6-spatial_analysis/wc2.1_30s_elev.zip")

