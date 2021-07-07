#'@title:  "Nulls networks altitude"
#'@author: "R Alcala"
#'@date:   "07/06/2021"
#'@output: "Nulls"
#'
#-- loading functions 
# p value calculator from z-scores
convert.z.score <- function(val){
  y=nrow(val)
  x = round(unlist(as.list(val[,-1])), 4)
  x = -abs(x)
  for (i in 1:length(x)){
    x[[i]] <- 2*pnorm(x[i])
  }
  y <- matrix(x, y, 3)
  res = as.data.frame(cbind(val[,1],y))
  colnames(res) <- colnames(val)
  res
}
#'
unlink(".RData")
setwd("~/Dropbox (UFL)/Alcala_Briseno-Garrett/++Papa_virome/+papa/3-results/")
load("PPV-ALL3-bipartitenetworkJul06.RData")
#------------ RESULTS
dimdatb.net.table[[1]]
dimdatb.net.table[[2]]
dimdatb.net.table[[3]]

#--
library(picante)
reps=100
#list 
findswap = cindswap = hindswap = list()
##---------------------------------------
# frozen land
fdat <- vina.ml[[1]]
femp.nest = nested(fdat, method="wine")
femp.mod <- computeModules(fdat)
femp.mod@likelihood

## Null hypothesis for independent swap:  
for(i in 1:reps){
  print(i)
  findswap[[i]] = randomizeMatrix(fdat, null.model="independentswap", iterations = 100)
}

fmod.null = fnest.null = fconn.null = list()

for(i in 1:reps){
  print(i)
  fnest.null[[i]] <- nested(findswap[[i]], method = "wine")
  fmod.null[[i]] <- computeModules(findswap[[i]])
}

fned.l = fmod.l = list()
for(i in 1:reps){
  print(i)
  fmod.l[[i]]  <- fmod.null[[i]]@likelihood
  fned.l[[i]] <- fnest.null[[i]]
}

fmod.null  <- (femp.mod@likelihood - mean(unlist(fmod.l))) / sd(unlist(fmod.l))
fnest.null  <- (femp.nest - mean(unlist(fned.l))) / sd(unlist(fned.l))

fval <- cbind(empirical=c(femp.mod@likelihood, femp.nest), null = c(fmod.null, fnest.null))
(fronzen.pval <- convert.z.score(fval)[1:2])

# plotting modularity 
plotModuleWeb(femp.mod, weighted = TRUE, plotModules = TRUE, labsize = 0.5)

#------- cold land
cdat <- vina.ml[[2]]

cemp.nest = nested(cdat, method="wine")
cemp.mod <- computeModules(cdat)
cemp.mod@likelihood
#--
# Null hypothesis for independent swap:  
for(i in 1:reps){
  print(i)
  cindswap[[i]] = randomizeMatrix(cdat, null.model="independentswap", iterations = 100)
}

cmod.null = cnest.null = cconn.null = list()

for(i in 1:reps){
  print(i)
  cnest.null[[i]] <- nested(cindswap[[i]], method = "wine")
  cmod.null[[i]] <- computeModules(cindswap[[i]])
}

cned.l = cmod.l = list()
for(i in 1:reps){
  print(i)
  cmod.l[[i]]  <- cmod.null[[i]]@likelihood
  cned.l[[i]] <- cnest.null[[i]]
}

cmod.null  <- (cemp.mod@likelihood - mean(unlist(cmod.l))) / sd(unlist(cmod.l))
cnest.null  <- (cemp.nest - mean(unlist(cned.l))) / sd(unlist(cned.l))

cval <- cbind(empirical=c(cemp.mod@likelihood, cemp.nest), null = c(cmod.null, cnest.null))
(coldland.pval <- convert.z.score(cval)[1:2])

# plotting modularity 
plotModuleWeb(cemp.mod, weighted = TRUE, plotModules = TRUE, labsize = 0.5)

#------------- hot land
hdat <- vina.ml[[3]]

hemp.nest = nested(hdat, method="wine")
hemp.mod <- computeModules(hdat)
hemp.mod@likelihood
#--
# Null hypothesis for independent swap:  
hindswap = list()
for(i in 1:reps){
  print(i)
  hindswap[[i]] = randomizeMatrix(hdat, null.model="independentswap", iterations = 100)
}

hmod.null = hnest.null = hconn.null =list()

for(i in 1:reps){
  print(i)
  hnest.null[[i]] <- nested(hindswap[[i]], method = "wine")
  hmod.null[[i]] <- computeModules(hindswap[[i]])
}

hned.l = hmod.l = list()
for(i in 1:reps){
  print(i)
  hmod.l[[i]]  <- hmod.null[[i]]@likelihood
  hned.l[[i]] <- hnest.null[[i]]
}

hmod.null  <- (hemp.mod@likelihood - mean(unlist(hmod.l))) / sd(unlist(hmod.l))
hnest.null  <- (hemp.nest - mean(unlist(hned.l))) / sd(unlist(hned.l))

hval <- cbind(empirical=c(hemp.mod@likelihood, hemp.nest), null = c(hmod.null, hnest.null))
(hotland.pval <- convert.z.score(hval)[1:2])

# plotting modularity 
plotModuleWeb(hemp.mod, weighted = TRUE, plotModules = TRUE, labsize = 0.5)


pval.res <- rbind(fronzen.pval,
      coldland.pval,
      hotland.pval)
#------- writing files
z=1 # 2 or 3
zones[z]
write.table(dimdatb.net.table[[z]], paste0("PPV_",zones[z], "_network-metrics_", format(Sys.time(), "%b%d"), ".tbl"), sep = ",")


#--
#
#


tic("Node metrics")
print("Processing node metrics")
datb.sp.table  <- specieslevel(vina.ml[[1]])
nodelev1 <- as.data.frame(datb.sp.table$`higher level`)
nodelev2 <- as.data.frame(datb.sp.table$`lower level`)
write.table(nodelev1, paste0("papa", zones[z],format(Sys.time(), "%b%d"),"_node-metrics_hl.tbl"))
write.table(nodelev2, paste0("papa", zones[z],format(Sys.time(), "%b%d"),"_node-metrics_ll.tbl"))
toc()

#----------------
#------- graph
dat.mat <- dat.mat[which(rownames(vina.ml[[1]]) %in% rownames(datb.sp.table$`lower level`)),]
dat.mat[is.na(dat.mat)] <- 0
# write.table(dat.mat, paste0("papa","zones[z]", format(Sys.time(), "%b%d"),"_data_matrix.tbl"))
tic("Plotting graph")
print("Creating graph")
virome <-  graph.incidence(na.omit(dat.mat), weighted=T)
#------- Network
# Prepraring colors
a <- rowSums(dat.mat)
snode <- data.frame(Species=names(a), Coverage=as.integer(a)) #normalized RPKMs
snode <- snode[snode$Coverage > 0,]
snodes = rbind(
          snode, Species=rep("nodes", ), Coverage=rep(100, n)))
# aading attributes
V(virome)$type
V(virome)$name <- c(V(virome)$name[1:length(V(virome)$type[V(virome)$type == "FALSE"])],
                    rep("", length(V(virome)$type[V(virome)$type == "TRUE"])))

V(virome)$color <-  c(rep("#FFFB95", length(V(virome)$type[V(virome)$type == "FALSE"])), rep("#92C5FC", length(V(virome)$type[V(virome)$type == "TRUE"])))
V(virome)$xx <- log2(snode$Coverage) # coverage size
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




# Saving R image
print("Saving R image")
save.image(paste0("PPV-ALL3", "-bipartitenetwork",format(Sys.time(), "%b%d"),".RData"))
toc()
