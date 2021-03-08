#'@title:  "rpkms-x-vsc loader"
#'@author: "R Alcala"
#'@date:   "15/02/2021"
#'@output: "rpkmsx-vsc merged"
#'
unlink(".RData")

library(plyr)

#--------
setwd("/Users/ricardoi/Dropbox (UFL)/Alcala_Briseno-Garrett/+++Sweetpotato_virome/+Sweetpotato_virome/rpkm_resultsFeb2/")


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

# add name
sweetpotato_virome <- rbind.fill(ls) 


write.csv(sweetpotato_virome, "sweetpotato_virome_vsc-rpkmx_ViNAtq_Feb24.csv")

ppv= peruvian_potato_virome
# how to read all csv's 
# dat.x <- ddply(ppv, .(IDs, Species), summarise, RPKM=mean(RPKM_mean))
# head(dat.x)

