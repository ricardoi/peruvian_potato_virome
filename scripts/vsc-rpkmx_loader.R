setwd("/Users/ricardoi/Dropbox (UFL)/Alcala_Briseno-Garrett/++Papa_virome/+papa/3-results/vsc-rpkmx")


files = list.files(pattern="*.csv")
files
#
ls = list(NULL)
for (i in seq_along(files)){
  x <- read.csv(files[i], as.is=T)
  x$IDs <- rep(strsplit(files[i], "_")[[1]][1], nrow(x))
  ls[[i]] = x
}
ls[[1]]

