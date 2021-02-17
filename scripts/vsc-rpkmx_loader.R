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

peruvian_potato_virome <- rbind.fill(ls) 


write.csv(peruvian_potato_virome, "peruvian_potato_virome_vsc-rpkmx_ViNAtq_Feb11.csv")

ppv= peruvian_potato_virome
# how to read all csv's 
# dat.x <- ddply(ppv, .(IDs, Species), summarise, RPKM=mean(RPKM_mean))
# head(dat.x)

