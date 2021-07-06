#'@title:  "Bipartite networks altitude"
#'@author: "R Alcala"
#'@date:   "13/03/2021"
#'@output: "Metadata"
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


#----- loading arguments 
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
# ggplot(datc, aes(Genus))+
#   geom_bar()

counts <- datc[which(datc$Genus == "Torradovirus"),]
unique(counts$pSpecies)

alluvial(datc[,c(1,2,3,5)], freq=datc$RPKM_mean,
         #hide = datc$length == 0,
         col = datc$cols,
         border = datc$cols, 
         alpha = 0.7,
         blocks = FALSE,
         ordering = list(NULL, NULL,NULL, order(as.factor(datc$Realm))), 
         # change NULL to order them
         cex =0.8
)

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
z=2 # 1 to 2 or 3
zones[z]
#--------------------------------------
# make an elseif statment to choose cluster number if specified, or full otherwise 
print(paste("Peruvian potato - altitudinal regions ::", zones[z]))
#----- 
k1vina <- virome2[which(virome2$altzones == zones[z]),]
# !comenta  for full network 
# k1vina = kvina # This line is to compute the full network ## Comment always 
k1meta <- meta[which(meta$altzones == zones[z]),]
# !comenta  for full network 
# k1meta = meta # This line is to compute the full network ## Comment always 
k1meta <- k1meta[which(k1meta$SampleID %in% unique(k1vina$IDs)),]
head(k1vina)

#subset >10
k1vina.10 <- k1vina[which(k1vina$Bases_mean > 0),]

pdf(paste0("PPV_", zones[z], "_RPKM_Genus_boxplot",format(Sys.time(), "%b%d"), ".pdf"),
    width = 18, # The width of the plot in inches
    height = 10) # The height of the plot in inches
ggplot(data = k1vina.10, aes(x = Genus, y = log2(Coverage_mean)))+
  geom_boxplot()+
  # stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")+
  geom_jitter(shape=1, aes(colour = Realm), alpha = 0.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#--
k1vina.att <- ddply(k1vina.10, .(Acronym, IDs), summarise, Coverage=mean(RPKM_mean))
k1vina.m <- tidyr::spread(k1vina.att, IDs, Coverage,  drop=TRUE , fill = 0)
rownames(k1vina.m) <- k1vina.m$Acronym
datm <- k1vina.m[-c(1)]
n=length(datm)
paste("MESSAGE:: The lenght of the sample locations is ", length(datm), " and the data matrix is ", n,
      ", then dimesions are equal? The answer is ", length(unique(k1vina$IDs)) == length(datm), sep = "")

#--- Normalization
print(paste("normalizing dataset", zones[z]))
dat.mat <- datm#/colSums(datm)*100
dat.mat = round(dat.mat, digits = 0)
dim(dat.mat)
dat.mat

dat.mat.log <- log2(dat.mat)
dat.mat.log[dat.mat.log==-Inf] <- 0
boxplot(dat.mat.log)


#------- Bipartite networks -------
print("Bipartite Network Analsysis began ... this may take a while")
# dat.mat = round(datm, digits = 0)
pdf(paste0(zones[z],"_bipartitenetwork_",format(Sys.time(), "%b%d"), ".pdf"),
    width = 15, # The width of the plot in inches
    height = 15) # The height of the plot in inches
# plot
plotweb(sortweb(dat.mat.log, sort.order="inc"), method="cca", abuns.type="additional",
        col.interaction="gray", text.low.col="gray1", text.rot=90)   

dev.off()
#pdf(paste0(k,"-kcluster_matrixnetwork",format(Sys.time(), "%b%d"), ".pdf"),
#    width = 15, # The width of the plot in inches
#    height = 15) # The height of the plot in inches
visweb(sortweb(dat.mat,sort.order="inc"), type="diagonal", labsize=3,
       square="interaction", text="none", textsize = 4,circles=FALSE, frame=FALSE)# 
#dev.off()
# dat.mat0 = ifelse(dat.mat > 1, 0, 1)

dat.df <- data.frame( Species = rownames(as.data.frame(rowSums(dat.mat))),
                      count =  as.numeric(rowSums(dat.mat)))
dat.bin = dat.mat
dat.bin[dat.bin >= 1 ] <- 1
dat.freq <- data.frame(Species = rownames(as.data.frame(rowSums(dat.bin))),
                       count =  as.numeric(rowSums(dat.bin)))

# dat.freq=dat.freq[dat.freq$count > 25,]

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

dat.df = dat.df[dat.df$count > 0,]
dat.mat <- dat.mat[rownames(dat.mat) %in% dat.df$Species,]
dim(dat.mat)


# Network information
tic("Network metrics")
print("Processing network metrics")
dimdatb.net.table <- networklevel(dat.mat)
write.table(dimdatb.net.table, paste0("papa-", zones[z],"_network-metrics.tbl", format(Sys.time(), "%b%d")), sep = ",")
toc()
tic("Node metrics")
print("Processing node metrics")
datb.sp.table  <- specieslevel(dat.mat)
nodelev1 <- as.data.frame(datb.sp.table$`higher level`)
nodelev2 <- as.data.frame(datb.sp.table$`lower level`)
write.table(nodelev1, paste0("papa", zones[z],format(Sys.time(), "%b%d"),"_node-metrics_hl.tbl"))
write.table(nodelev2, paste0("papa", zones[z],format(Sys.time(), "%b%d"),"_node-metrics_ll.tbl"))
toc()
# tic("Nestedness")
# datb.nest.table <- nestedcontribution(dat.mat)
# toc()
#------- graph
dat.mat <- dat.mat[which(rownames(dat.mat) %in% rownames(datb.sp.table$`lower level`)),]
write.table(dat.mat, paste0("papa","zones[z]", format(Sys.time(), "%b%d"),"_data_matrix.tbl"))
tic("Plotting graph")
print("Creating graph")
virome <-  graph.incidence(dat.mat, weighted=T)
#------- Network
# Prepraring colors
a <- rowSums(dat.mat)
snode <- data.frame(Species=names(a), Coverage=as.integer(a)) #normalized RPKMs
snodes = rbind(snode, data.frame(Species=rep("nodes", n), Coverage=rep(100, n)))
# aading attributes
V(virome)$type
V(virome)$name <- c(V(virome)$name[1:length(V(virome)$type[V(virome)$type == "FALSE"])],
                    rep("", length(V(virome)$type[V(virome)$type == "TRUE"])))

V(virome)$color <-  c(rep("#FFFB95", length(V(virome)$type[V(virome)$type == "FALSE"])), rep("#92C5FC", length(V(virome)$type[V(virome)$type == "TRUE"])))
V(virome)$xx <- log(snodes$Coverage) # coverage size
V(virome)$width <- c((datb.sp.table$`lower level`$species.strength)+6, (dim(dat.mat)[2]/(datb.sp.table$`higher level`$degree)*2))
show_col(unique(V(virome)$color))
# 
submeta <- meta[which(meta$SampleID %in% colnames(dat.mat)),]
legend <- unique(cbind(ColorCode=submeta$colors, Region=as.character(submeta$Country)))
# Edges names and attributes
rbPal <- colorRampPalette(c("#f5f5f5", "#483C33")) 
# rbPal <- colorRampPalette(c("yellow", "red3")) 

E(virome)$width <- (E(virome)$weight)
E(virome)$width[E(virome)$width <= 1] <- NA
E(virome)$width <- log2(E(virome)$width)
E(virome)$width[E(virome)$width == 0] <- NA
E(virome)$width <- (E(virome)$width)/max(E(virome)$width[!is.na(E(virome)$width)])
E(virome)$color <- rbPal(4)[as.numeric(cut(E(virome)$width, breaks = 4))]
show_col(unique(E(virome)$color))


# E(virome)$color <- rbPal(10)[as.numeric(cut(E(virome)$width, breaks = 10))]
shapes = c(rep("circle", length(V(virome)$type[V(virome)$type == "FALSE"])), rep("square", length(V(virome)$type[V(virome)$type == "TRUE"])))
# # Network
pdf(paste0(zones[z],"_bipartitenetwork-mds_",format(Sys.time(), "%b%d"), ".pdf"),
    width = 15, # The width of the plot in inches
    height = 15) # The height of the plot in inches
plot(virome, edge.arrow.size=1, vertex.shape=shapes, vertex.size=V(virome)$width , 
          vertex.label.cex=1, vertex.label.color='black', vertex.frame.color="gray", 
     vertex.frame.color="gold",   edge.curved=F,  layout=layout_with_mds(virome)) 

# dim = 3, maxiter = vcount(virome)*10)) # maxiter=500 ,fineiter = 500, cool.fact=0.80, weight.edge.crossings = 1 - sqrt(edge_density(virome))))
# legend(x=-1.3, y=1.1, legend[,2], pch=21,
#        col=legend[,1], pt.bg=legend[,1], pt.cex=1, cex=.8, bty="n", ncol=1)
dev.off()
toc()
#----------------------- TEST FOR RAREFACTION/EXTRAPOLATION
#---
# data(BCI) 
# print("Rarefaction/Extrapolation")
# dim(dat.mat)
# dat.mat[1:5,1:5]
# S <- specnumber(t(dat.mat)) # observed number of species
# (raremax <- min(rowSums(t(dat.mat))))
# Srare <- rarefy(t(dat.mat), raremax)
# plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
# abline(0, 1)
# pdf(paste0(k,"-kcluster_rarefaction-vegan_",format(Sys.time(), "%b%d"), ".pdf"),
#    width = 15, # The width of the plot in inches
#    height = 15) # The height of the plot in inches
# rarecurve(t(dat.mat), step = 20, sample = raremax, col = "blue", cex = 0.6)
#dev.off()
#---
## Truncating values to a max of 1000 for visualization 
# con <- ifelse(dat.mat <= 1000, print(1), Inf)
# dat.mat2=as.matrix(dat.mat) * as.matrix(con)
# dat.mat2[sapply(dat.mat2, is.infinite)] <- 1000
# dat.mat2[1:5, 1:5]

out2 <- iNEXT(t(dat.mat), q = 0, datatype = "abundance", se = T, nboot = 20)
pdf(paste0("papa", zones[z],"rarefaction-iNEXT_",format(Sys.time(), "%b%d"), ".pdf"),
    width = 15, # The width of the plot in inches
    height = 15) # The height of the plot in inches
ggiNEXT(out2, type = 1, color.var="none")+
  theme_bw()
dev.off()

# Saving R image
print("Saving R image")
save.image(paste0("papa", zones[z],"-bipartitenetwork",format(Sys.time(), "%b%d"),".RData"))
toc()
