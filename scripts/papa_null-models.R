#'@title:  "Null models - species associationº"
#'@author: "R Alcala"
#'@date:   "04/01/2021"
#'@output: "null models :: z-scores"
#'
unlink(".RData")
setwd(dir = ".")
#--------------------------------------------------------------------------------------
#---- Libraries ----
library(picante)
library(bipartite)
library(tictoc)
#---- functions ----
source("~/ViNA/bin/functions/func_nullmodels.R")	
print("libraries loaded and functions sourced")
#--------------------------------------------------------------------------------------
# Null model analysis of species association abundance data
# The authors created a series of structured matrices by altering the random matrices 
# to incorporate patterns of pairwise species segregation and aggregation.
# To test random and structure matrices to determine which tests had low Type I error rates 
# and good power for detecting segregated and aggregated species distribution.
# bundance matrices, analyzed with an appropriate null model, may be a powerful tool for 
# quantifying patterns of species segregation and aggregation

# Ulrich, W., & Gotelli, N. J. (2010). Null model analysis of species associations 
# using abundance data. Ecology, 91(11), 3384–3397.
# Note: The "matrix" software used was developed in Fortran for Windows, here I am using 
#       the package "picante" to randomized the matrices.
###############################################################################################
################################################################################################
#################################################################################################
# #------------------------------------------------------------------------------
# #-------- Functions
# # Bipartite Network metric Null calculator
bine.trics <- function(net){
  conn.null = nestw.null =ISA.null= SA.null = H2.null = list()
  for(i in 1:reps){
    conn.null[[i]] <- net[[i]][1]
    #nesth.null[[i]] <- net[[i]][8]
    nestw.null[[i]] <- net[[i]][8]
    ISA.null[[i]] <- net[[i]][10]
    SA.null[[i]] <- net[[i]][11]
    H2.null[[i]] <- net[[i]][18]
  }
  list(conn.null, nestw.null, ISA.null, SA.null, H2.null)
}

#--------------------------------------------------------------------------------------
# Z scores calculator
binet.zsc <- function(dat, null){
  conn.null.test  <- (dat[1] -mean(unlist(null[[1]]))) / sd(unlist(null[[1]]))
  #nesth.null.test<- (dat[7] -mean(unlist(null[[2]]))) / sd(unlist(null[[2]]))
  nestw.null.test <- (dat[8] -mean(unlist(null[[2]]))) / sd(unlist(null[[2]]))
  ISA.null.test   <- (dat[10] -mean(unlist(null[[3]]))) / sd(unlist(null[[3]]))
  SA.null.test    <- (dat[11] -mean(unlist(null[[4]]))) / sd(unlist(null[[4]]))
  H2.null.test    <- (dat[18]-mean(unlist(null[[5]]))) / sd(unlist(null[[5]]))
  c(conn.null.test, nestw.null.test, ISA.null.test, SA.null.test, H2.null.test)
}
#--------------------------------------------------------------------------------------
# p values calculator from z-score with stars ***= >0.001, **= >0.01, *= >0.5, ns < 0.5, na = not aplicable
pvalbystars.cal <- function(x){
  y <- matrix(0, nrow(x), ncol(x))
  for( i in 1:nrow(x)){
    for( j in 2:ncol(x)){
      if ( is.na(x[i,j] )) {
        y[i,j] <- "NA"
      } else if ( x[i,j] < -3.291 | x[i,j] > 3.291 ) {
        y[i,j] <- "***"
      } else if ( x[i,j] < -2.567 | x[i,j] > 2.567 ) {
        y[i,j] <- "**"
      } else if (x[i,j] < -1.96 | x[i,j] > 1.96 ) {
        y[i,j] <- "*"
      } else if (x[i,j] > -1.96 | x[i,j] < 1.96) {
        y[i,j] <- "ns"
      } else {
        print("NA")
      }
    }
  }
  res = as.data.frame(cbind(x[,1],y[,2:ncol(x)]))
  colnames(res) <- colnames(x)
  res
}
#--------------------------------------------------------------------------------------
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
################################################################################################
################################################################################################
################################################################################################
# #--------------------------------------------------------------------------------------
#------ Methods 
# 100 x n x m reshufflings
#
# FILE <- commandArgs(TRUE)
FILE <- "papacoldland-bipartitenetworkMar15.RData"
print(FILE)
#-----
load(file = FILE)
dat=dat.mat
print(paste("Matrix loaded::", dim(dat)[1], "x", dim(dat)[2]))
dat[1:5,1:5]
reps=100
# empirical <- networklevel(dat.p)
empirical <- networklevel(dat, ISAmethod="Bluethgen",  SAmethod = "Bluethgen")
emp.nest = nested(dat, method="wine")
emp.nest
#---------------------------------------
#----- Null hypothesis for richness:  
rich = list()
for(i in 1:reps){
  print(i)
   rich[[i]] = randomizeMatrix(dat, null.model="richness", iterations = 100)
}
print(paste("richness randomization for", FILE, ":: done"))
#---------------------------------------
#----- Null hypothesis for frequency:  
freq = list()
for(i in 1:reps){
  print(i)
   freq[[i]] = randomizeMatrix(dat, null.model="frequency", iterations = 100)
}
print(paste("frequency randomization for", FILE, ":: done"))
#----------------------------------------
#----- Null hypothesis for independent swap:  
indswap = list()
for(i in 1:reps){
  print(i)
  indswap[[i]] = randomizeMatrix(dat, null.model="independentswap", iterations = 100)
}

#--------------------------------------------------------------------------------------
# # Null hypothesis for trial swap:  
# triswap = list()
# for(i in 1:reps){
#   print(i)
#   triswap[[i]] = randomizeMatrix(dat.p, null.model="trialswap", iterations = 100)
# }
#--------------------------------------------------------------------------------------
#----- Calculate a uniform distribution of the Null models at network level
null1.network =list()
for(i in 1:reps){
  print(i)
  null1.network[[i]] <- networklevel(rich[[i]], ISAmethod="Bluethgen",  SAmethod = "Bluethgen")
}
#-----
null2.network =list()
for(i in 1:reps){
  print(i)
  null2.network[[i]] <- networklevel(freq[[i]], ISAmethod="Bluethgen",  SAmethod = "Bluethgen")
}
#-----
null3.network =list()
for(i in 1:reps){
  print(i)
  null3.network[[i]] <- networklevel(indswap[[i]], ISAmethod="Bluethgen",  SAmethod = "Bluethgen")
}
names <- names(empirical); 
names
#----------------------------------------------------------------------
#------- Modularity 
res <- computeModules(dat)
plotModuleWeb(res)

#------- Null for richness
null1.module =list()
for(i in 1:reps){
  print(i)
  null1.module[[i]] <- computeModules(rich[[i]])
}
#--
null1.mod = list()
for (i in 1:reps){
  print(i)
  null1.mod[[i]] <- null1.module[[i]]@likelihood
}
mean(unlist(null1.mod))
plotModuleWeb(null1.module[[3]])

mod.null1.test  <- (res@likelihood - mean(unlist(null1.mod))) / sd(unlist(null1.mod))

#------- Null for frequency
null2.module =list()
for(i in 1:reps){
  print(i)
  null2.module[[i]] <- computeModules(freq[[i]])
}
#--
null2.mod = list()
for (i in 1:reps){
  print(i)
  null2.mod[[i]] <- null2.module[[i]]@likelihood
}
mean(unlist(null2.mod))
plotModuleWeb(null2.module[[1]])

mod.null2.test  <- (res@likelihood - mean(unlist(null2.mod))) / sd(unlist(null2.mod))

#------- Null for indswap
null3.module =list()
for(i in 1:reps){
  print(i)
  null3.module[[i]] <- computeModules(indswap[[i]])
}
#--
null3.mod = list()
for (i in 1:reps){
  print(i)
  null3.mod[[i]] <- null3.module[[i]]@likelihood
}
mean(unlist(null3.mod))
plotModuleWeb(null3.module[[1]])

mod.null3.test  <- (res@likelihood - mean(unlist(null3.mod))) / sd(unlist(null3.mod))
#--------------------------------------------------------------------------------------------------------
#----- P A R A L L E L 
# Calculate a uniform distribution of the Null models at network level
library(foreach)
library(doParallel)
#setup parallel backend to use many processors
cores=detectCores()
cores
cl <- makeCluster(8)#cores[1]-2) #not to overload your computer
registerDoParallel(cl)

null1.network <-  foreach(i = 1:reps, .packages= "bipartite", .combine = rbind) %dopar% {
                    networklevel(rich[[i]], ISAmethod="Bluethgen",  SAmethod = "Bluethgen")
}

null2.network <-  foreach(i = 1:reps, .packages= "bipartite", .combine = rbind) %dopar% {
                    networklevel(freq[[i]], ISAmethod="Bluethgen",  SAmethod = "Bluethgen")
}

null3.network <-  foreach(i = 1:reps, .packages= "bipartite", .combine = rbind) %dopar% {
                    networklevel(indswap[[i]], ISAmethod="Bluethgen",  SAmethod = "Bluethgen")
}
#--------------------------------------------------------------------------------------
# Calculate a uniform distribution of the Null models for modularity 

null1.module <- foreach(i = 1:reps, .packages = "bipartite", .combine=rbind) %dopar% {
          test <- computeModules(rich[[i]])
  test@likelihood
}

#--------------------------------------------------------------------------------------
#----- Results
#Bipartite network metric calculation
null.1 <- bine.trics(null1.network)
null.2 <- bine.trics(null2.network)
null.3 <- bine.trics(null3.network)

#----- Checking normal data:
#-------------------------------
#pdf(file= "null.qqplot.pdf")
par(mfrow=c(5,3))
for(i in 1:length(null.1)){
  qqplot(c(1:1000), unlist(null.1[[i]]))
  qqplot(c(1:1000), unlist(null.2[[i]]))
  qqplot(c(1:1000), unlist(null.3[[i]]))
}
# dev.off()
#---
pdf(file = "null.histall.pdf", width = 10, height = 18)
par(mfrow=c(6,3))
hist(unlist(null1.mod))
hist(unlist(null2.mod))
hist(unlist(null3.mod))
for(i in 1:length(null.1)){
  hist(unlist(null.1[[i]]))
  hist(unlist(null.2[[i]]))
  hist(unlist(null.3[[i]]))
}
dev.off()
#-------------------------------
# Result z-scores by 
res.null.1 <- binet.zsc(empirical, null.1)
res.null.2 <- binet.zsc(empirical, null.2)
res.null.3 <- binet.zsc(empirical, null.3)

val <- cbind(empirical=c(empirical[c(1, 8, 10, 11, 18)]), null1=res.null.1, null2=res.null.2, null3=res.null.3)
#-------------------------------
# Modularity
modularity.null <- cbind(res@likelihood, mod.null1.test, mod.null2.test, mod.null3.test)
values <- rbind(modularity.null, val)
rownames(values) <- c("modularity", rownames(values)[-1])
colnames(values) <- c("observerd", "Null.model.1", "Null.model.2", "Null.model.3")
#-------------------------------
# Conver z-scores to p values
z.score.results = convert.z.score(values)
z.score.results
#-------------------------------
# Conver z-scores to p values by stars
z.score.resultsbystars = pvalbystars.cal(round(values, 2))
z.score.resultsbystars


write.csv(z.score.results, "Nestedneess_null-model/results_networks.csv")
write.csv(z.score.resultsbystars, "Nestedneess_null-model/results_networksstars.csv")


# results.all
# results.cost
# results.depre

# z-score: 
# 1.96 = 0.5 pvalue
# 2.576 = 0.1 pvalue
# 3.291 = 0.001 pvalue
# round(pnorm(1.96) - pnorm(-1.96), 2) - 1
# Positive reltionship between Species strenght and species degree
sp.lev <- specieslevel(dat.mat)
sp.lev
#lineal
plot(sp.lev$`lower level`$species.strength, sp.lev$`lower level`$degree)
sp.lm <- lm(sp.lev$`lower level`$degree ~ sp.lev$`lower level`$species.strength)
abline(sp.lm)
summary(sp.lm)

# quadratic
plot(sp.lev$`lower level`$species.strength, sp.lev$`lower level`$degree)
sp.qd <- lm(sp.lev$`lower level`$degree ~ 1 + sp.lev$`lower level`$species.strength +
              I((sp.lev$`lower level`$species.strength)^2))
points(sp.lev$`lower level`$degree, predict(sp.qd), type="l")
summary(sp.qd)

#-----
print("Saving R image")
file <- strsplit(FILE, "-b")[[1]][1]
save.image(paste0(file,"null_models",format(Sys.time(), "%b%d"),".RData"))

