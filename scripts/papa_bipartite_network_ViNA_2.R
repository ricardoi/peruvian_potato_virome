#'@title:  "Bipartite networks altitude"
#'@author: "R Alcala"
#'@date:   "13/03/2021"
#'@output: "Boxplot RPKMs"
#'
unlink(".RData")
setwd("~/Dropbox (UFL)/Alcala_Briseno-Garrett/++Papa_virome/+papa/3-results/")
#---- Libraries ----
library(plyr)
library(dplyr)
library(tidyr)
library(igraph)
library(bipartite)
library(reshape2)
library(viridis)
library(scales)
library(ggplot2)
library(vegan)
library(iNEXT)
library(tictoc)
library(gridExtra)
library(openxlsx)
#----- loading functions

#---------- Setting working directory ---------- 
tic("Total time computed")
{
  tic("Loading and formating")
  print("Loading data and metadata")
  
  # ----- Adding metadata --------
  metadata <- read.xlsx("2_Viroma_web_180319_final.xlsx", sheet = 2)
  metadata <- metadata %>% select_if(~!all(is.na(.)))
  metadata
  
  unique(metadata$Localidad)
  
  length(metadata$`Altitude.(masl)`)
  hl <- metadata[which(metadata$`Altitude.(masl)` < 1000),] # hot land Hoehen index
  hls<- c(table(hl$Department), Total= sum(table(hl$Department)))
  cl <- metadata[which(metadata$`Altitude.(masl)` > 2000 & metadata$`Altitude.(masl)` <3400), ] # cold land Hoehen index
  cls <- c(table(cl$Department), Total= sum(table(cl$Department)))
  fl <- metadata[which(metadata$`Altitude.(masl)` > 3400),] # frozen land Hoehen index
  fls <- c(table(fl$Department), Total= sum(table(fl$Department)))
  
  metadata$altzones <- ifelse(metadata$`Altitude.(masl)` < 1000, "hotland",
                              ifelse(metadata$`Altitude.(masl)` > 1000 & metadata$`Altitude.(masl)` <2000, "templand", 
                                     ifelse(metadata$`Altitude.(masl)` > 2000 & metadata$`Altitude.(masl)` <3400, "coldland",
                                            ifelse(metadata$`Altitude.(masl)` > 3400 & metadata$`Altitude.(masl)` <5000, "frozenland",
                                                   ifelse(metadata$`Altitude.(masl)` > 5000, "snowline", "out of range")))))
  
  
  meta <-metadata[c("CIP_Code", "Department", "Province", "Localidad",
                    "Latitude", "Longitude", "Altitude.(masl)", "Tamaño.campo.(m2)",
                    "Cultivar/.especie","Origen.de.semilla", "altzones")]
  
  #---- Peruvian Potato Virome
  peruvian_potato_virome <- read.csv("peruvian_potato_virome_vsc-rpkmx_ViNAtq_ed_Mar15.csv", 
                                     stringsAsFactors = F)[-1]
  virome = merge(peruvian_potato_virome, meta[c("CIP_Code","Altitude.(masl)","altzones")], 
                 by.x="IDs", by.y="CIP_Code")
  
  virome[1:5,]
  dim(virome)
  
  # first plot
  plot(sort(metadata$`Altitude.(masl)`))
  abline(h =c(1000, 2000, 3400), col = "red")
  
  #---
  #---- Subset a dataframes ---- 
  virome$Family <- gsub("^ ", "", virome$Family)
  virome$Genus <- gsub("^ ", "", virome$Genus)
  virome$Species <- gsub("^ ", "", virome$Species)
  virome$pSpecies <- gsub("^ ", "", virome$pSpecies)
  virome$Acronym <- gsub("^ ", "", virome$Acronym)
  as_tibble(virome)
  
  table(virome$Family)
  table(virome$Genus)
  table(virome$Species)
  # Removing no plant virus hits
  cond <- virome$Genus %in% c("Unassigned virus", "Trichovirus", "Mitovirus", " Hypovirus","Pospiviroid", " NA")
  virome <- virome[!cond,]
  # summary(virome$Length_mean); summary(virome$RPKM_mean)
  condNA <- as.numeric(virome$Length_mean)
  condNA[is.na(condNA)] <- 0
  virome$Length_mean <- condNA
  
  
  # subsetting for plotting lenth and RPKMs 
  cond1= as.numeric(virome$Length_mean) > 50
  cond2= as.numeric(virome$Length_mean) < 50
  
  virome2 <- virome[cond1,]
  kvina2 <- virome[cond2,]
}
toc()
#-------- Taxonomic plot
tic()
library(alluvial)
library(tidyverse)
library(reshape2)
# subsetting table
datc <- virome2 %>%
  select(Realm, Family, Genus, pSpecies, Acronym, altzones, RPKM_mean)
head(datc)

# adding colors
library(wesanderson)

RbPal <- viridis(c(length(unique(datc$Family))+40), option = "inferno")
# RbPal <- wes_palette("Zissou1", 100, type = "continuous")
datc <- datc %>%
  mutate( ss = paste(Family, Genus),
          cols = RbPal[ match(ss, sort(unique(ss))) ] 
  )
#----- Taxa Alluvla Plot 

plot(sort(table(datc$Family), decreasing = T), col="blue",las =2, cex = 2, lwd = 6)
plot(sort(table(datc$Genus), decreasing = T), col="blue",las =2, cex = 2, lwd = 6)
# plot(sort(table(datc$Acronym), decreasing = T), col="blue",las =2, cex = .5, lwd = 4)
# 

counts <- datc[which(datc$Genus == "Torradovirus"),]
unique(counts$pSpecies)
#)

#---------------
# Plot of mean contig length after removal of 50 nt
# pdf(paste0("papa","-nornRPKM-length",format(Sys.time(), "%b%d"), ".pdf"),
#     width = 15, # The width of the plot in inches
#     height = 15) # The height of the plot in inches
par(mfrow=c(1,2))
plot(virome$Length_mean, log2(virome$RPKM_mean),  col="black",
     xlab = "contig length", ylab = "log(mean RPKM)")
title("Normalized mean \n RPKM and length\n
      (thr = all)")
plot(virome2$Length_mean, log2(virome2$RPKM_mean),  col="black", pch = 3,
     xlab = "contig length", ylab = "log(mean RPKM)")
points(kvina2$Length_mean, log2(kvina2$RPKM_mean), col="red", pch = 1)
title("Normalized mean \n RPKM and length \n
      (thr > 50 nt)")
# dev.off()

toc()


zones <- unique(virome$altzones)
# z=1# 1 to 2 or 3
zones[]
#--------------------------------------
# make an elseif statment to choose cluster number if specified, or full otherwise 
print(paste("Peruvian potato - altitudinal regions ::", zones[]))
#----- 
f.vina <- virome2[which(virome2$altzones == zones[1]),]
c.vina <- virome2[which(virome2$altzones == zones[2]),]
h.vina <- virome2[which(virome2$altzones == zones[3]),]

f.meta <- meta[which(meta$altzones == zones[1]),]
c.meta <- meta[which(meta$altzones == zones[2]),]
h.meta <- meta[which(meta$altzones == zones[3]),]

# k1meta <- k1meta[which(k1meta$SampleID %in% unique(k1vina$IDs)),]

f.vmeta <- f.meta[which(f.meta$CIP_Code %in% unique(f.vina$IDs)),]
c.vmeta <- c.meta[which(c.meta$CIP_Code %in% unique(c.vina$IDs)),]
h.vmeta <- h.meta[which(h.meta$CIP_Code %in% unique(h.vina$IDs)),]

head(f.vmeta);head(c.vmeta);head(h.vmeta)

#subset >10
# k1vina.10 <- k1vina[which(k1vina$Bases_mean > 0),]
vina.10 <- list(
f.vina.10 <- f.vina[which(f.vina$Bases_mean > 0),],
c.vina.10 <- c.vina[which(c.vina$Bases_mean > 0),],
h.vina.10 <- h.vina[which(h.vina$Bases_mean > 0),]
)

#-- boxplot

f <- ggplot(data = f.vina.10, aes(x = IDs, y = log2(RPKM_mean)))+
  geom_boxplot()+
  # stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")+
  geom_jitter(shape=1, aes(colour = Realm), alpha = 0.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
c <- ggplot(data = c.vina.10, aes(x = IDs, y = log2(RPKM_mean)))+
  geom_boxplot()+
  # stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")+
  geom_jitter(shape=1, aes(colour = Realm), alpha = 0.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
h <- ggplot(data = h.vina.10, aes(x = IDs, y = log2(RPKM_mean)))+
  geom_boxplot()+
  # stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")+
  geom_jitter(shape=1, aes(colour = Realm), alpha = 0.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0("PPV_", "Genus_3zones", "_RPKM_IDs_boxplot",format(Sys.time(), "%b%d"), ".pdf"),
    width = 30, # The width of the plot in inches
    height = 14) # The height of the plot in inches

grid.arrange(f, c, h, ncol=1)

dev.off()


#-- Creating matrix

vina.ml = vina.m = vina.m10 = list()

for (i in seq_along(vina.10)){
  vina.m10[[i]] <- ddply(vina.10[[i]], .(Acronym, IDs), summarise, Coverage=mean(RPKM_mean))
}

for (i in seq_along(vina.m10)){
    mat <- tidyr::spread(vina.m10[[i]], IDs, Coverage,  drop=TRUE , fill = 0)
      rownames(mat) <- mat$Acronym
  vina.m[[i]] <- mat[-c(1)]
}


# n=length(datm)
# paste("MESSAGE:: The lenght of the sample locations is ", length(datm), " and the data matrix is ", n,
      # ", then dimesions are equal? The answer is ", length(unique(k1vina$IDs)) == length(datm), sep = "")

#--- Normalization
for (z in seq_along(vina.m)){
  print(paste("normalizing dataset", zones[z]))
    dat.mat <- vina.m[[z]]#/colSums(datm)*100
    dat.mat = round(dat.mat, digits = 0)
    print(dat.mat[1:5,1:5])
    # log 2 transformation
    dat.mat.log <- log2(dat.mat)
    dat.mat.log[dat.mat.log==-Inf] <- 0
  vina.ml[[z]] <- round(dat.mat.log, digits = 0)
}
#------- Bipartite networks -------
print("Bipartite Network Analsysis began ... this may take a while")
z=1 # 2 or 3 
#
pdf(paste0(zones[z],"_bipartitenetwork_",format(Sys.time(), "%b%d"), ".pdf"),
    width = 15, # The width of the plot in inches
    height = 15) # The height of the plot in inches
# plot
plotweb(sortweb(vina.ml[[z]], sort.order = "inc"), method = "cca", abuns.type= "additional",
        col.interaction = "gray", text.low.col = "gray1", text.rot = 90)   
dev.off()

#pdf(paste0(k,"-kcluster_matrixnetwork",format(Sys.time(), "%b%d"), ".pdf"),
#    width = 15, # The width of the plot in inches
#    height = 15) # The height of the plot in inches
visweb(sortweb(vina.ml[[z]], sort.order = "inc"), type="diagonal", labsize=3,
       square="interaction", text="none", textsize = 4,circles=FALSE, frame=FALSE)# 
#dev.off()
# dat.mat0 = ifelse(dat.mat > 1, 0, 1)

dat.df <- data.frame( Species = rownames(as.data.frame(rowSums(dat.mat))),
                      count =  as.numeric(rowSums(dat.mat)))
dat.bin = dat.mat
dat.bin[dat.bin >= 1 ] <- 1
dat.freq <- data.frame(Species = rownames(as.data.frame(rowSums(dat.bin))),
                       count =  as.numeric(rowSums(dat.bin)))

#- plotting  RPKMs and frequencies
plot1 <- ggplot(dat.df, aes(x = Species, y = sort(count)))+
  geom_bar(stat = "identity")+
  coord_flip()+
  geom_hline(yintercept=50, linetype="dashed", color = "red")+
  theme_minimal()
plot2 <- ggplot(dat.freq, aes(x = Species, y = sort(count)))+
  geom_bar(stat = "identity")+
  coord_flip()+
  geom_hline(yintercept=25, linetype="dashed", color = "red")+
  theme_minimal()

pdf(paste0(zones[z],"_incidence_w+b", format(Sys.time(), "%b%d"), ".pdf"),
    width = 30, # The width of the plot in inches
    height = 15) # The height of the plot in inches
grid.arrange(plot1, plot2, ncol=2)
dev.off()


#---
# dat.df = dat.df[dat.df$count > 0,] # Removing with 0's
# dat.mat <- dat.mat[rownames(dat.mat) %in% dat.df$Species,]
dim(dat.mat)

#------ Calculating networks with bipartite
library(doParallel)
registerDoParallel(cores= seq_along(vina.ml) )  
# Network information
tic("Network metrics")
print("Processing network metrics")

dimdatb.net.table <- foreach(W=seq_along(vina.ml)) %dopar% 
  networklevel(vina.ml[[W]])

