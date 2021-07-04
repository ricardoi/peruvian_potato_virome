#'@title:  "Virome Networks Metric"
#'@author: "R Alcala"
#'@date:   "06/24/2021"
#'@output: "Plots metrics vs regions"
#---------------------------- virome network metrics  --------------------------- 
setwd("~/git_local/sweetpotato_virome/data/netork_metrics/")

#----------------- Install & load libraries ------------------
library(tidyverse)
#library(dplyr)
#library(plyr)
# library(igraph )
#### NOTE

#------- Loading files
files = list.files(pattern="*network-metrics.csv")
#
ls=list()
for (i in 1:length(files)){
  x <- read.csv(files[i], as.is=T)
  # x$IDs <- rep(names$names[i], nrow(x))
  ls[[i]] =  x
}
ls[[1]]

nmts <- cbind(ls[[1]], ls[[2]], ls[[3]], ls[[4]], ls[[5]], ls[[6]], ls[[7]])
colnames(nmts) <- c(1:7)
nmts <- as.data.frame(t(nmts))
#####
names(nmts)
# Figure
# replace regions with no. of individuals and the dots indicate what region is
par(mfrow=c(3,4), mar=c(4,4,4,4))
plot(nmts$number.of.species.LL ~ nmts$number.of.species.HL, xlab = "No. of species (host)", ylab = "No. of species (virus)")
text(nmts$number.of.species.LL ~ nmts$number.of.species.HL, cex=0.9, font=2, pos=2)
     # main="host-virus species")
plot(nmts$nestedness ~ nmts$number.of.species.HL , ylab = "Nestedness", xlab = "No. of species (host)")
text(nmts$nestedness ~ nmts$number.of.species.HL, cex=0.9, font=2, pos=2)
# main="nestedness")
plot(nmts$`weighted nestedness` ~ nmts$number.of.species.HL, ylab = "weighted nestedness", xlab = "No. of species (host)") 
text(nmts$`weighted nestedness` ~ nmts$number.of.species.HL, cex=0.9, font=2, pos=2)
# main="weighted nestedness")
plot(nmts$NODF ~ nmts$number.of.species.HL, ylab = "NODF", xlab="No. of species (host)")
text(nmts$NODF ~ nmts$number.of.species.HL, cex=0.9, font=2, pos=2)
# main="NODF")
plot(nmts$`connectance` ~ nmts$number.of.species.HL, ylab = "connectance", xlab="No. of species (host)")
text(nmts$`connectance` ~ nmts$number.of.species.HL, cex=0.9, font=2, pos=2)
# main="connectance")
plot(nmts$`weighted connectance` ~ nmts$number.of.species.HL, ylab = "weighted connectance", xlab="No. of species (host)")
text(nmts$`weighted connectance` ~ nmts$number.of.species.HL, cex=0.9, font=2, pos=2)
# main="weighted connectance")
plot(nmts$`Fisher alpha` ~ nmts$number.of.species.HL, ylab = "Fisher alpha", xlab="No. of species (host)")
text(nmts$`Fisher alpha` ~ nmts$number.of.species.HL, cex=0.9, font=2, pos=2)
# main="Fisher alpha")
plot(nmts$`Shannon diversity` ~ nmts$number.of.species.HL, ylab = "Shannon diversity", xlab="No. of species (host)")
text(nmts$`Shannon diversity` ~ nmts$number.of.species.HL, cex=0.9, font=2, pos=2)
# main="Shannon diversity")
plot(nmts$`interaction evenness` ~ nmts$number.of.species.HL, ylab = "interaction evenness", xlab="No. of species (host)")
text(nmts$`interaction evenness` ~ nmts$number.of.species.HL, cex=0.9, font=2, pos=2)
# main="interaction evenness")
plot(nmts$`Alatalo interaction evenness` ~ nmts$number.of.species.HL, ylab = "Alatalo interaction evenness", xlab="No. of species (host)")
text(nmts$`Alatalo interaction evenness` ~ nmts$number.of.species.HL, cex=0.9, font=2, pos=2)
# main="Alatalo interaction evenness")
plot(nmts$C.score.HL ~ nmts$number.of.species.HL, ylab = "interaction evenness", xlab="No. of species (host)")
text(nmts$C.score.HL ~ nmts$number.of.species.HL, cex=0.9, font=2, pos=2)
# main="C-score HL")
plot(nmts$V.ratio.HL ~ nmts$number.of.species.HL, ylab = "interaction evenness", xlab="No. of species (host)")
text(nmts$V.ratio.HL ~ nmts$number.of.species.HL, cex=0.9, font=2, pos=2)
# main="V-score HL")

# dev.off()

#------ Node metrics HL
files2 = list.files(pattern="*metrics_hl.csv")
names <- strsplit(files2, split = "_")
#
for (i in 1:length(files2)){
  x <- read.csv(files2[i], as.is=T)
  x$IDs <- rep(names[[i]][1], nrow(x))
  ls[[i]] = x
  colnames(ls[[i]]) <- colnames(ls[[1]])
}

nohl <- rbind(ls[[1]], ls[[2]], ls[[3]], ls[[4]], ls[[5]], ls[[6]], ls[[7]])
class(nohl)
as_tibble(nohl)

ggplot(nohl, aes(IDs, degree))+ 
  geom_count(shape = 1, show.legend=T) +
  labs(subtitle="Higher level SSA-SPV", 
       y="Node degree", 
       x="Sweetpotato region", 
       title="Count Points")
#-
library(ggthemes)
# plot
#----- degree
ggplot(nohl, aes(IDs, degree))+
  geom_tufteboxplot(voffset = 0.0, hoffset = 0, median.type = "line") + 
  # geom_point()+
  stat_summary(fun = "mean", geom = "point", col = "black", size = 3, shape = 1)+
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  # theme_tufte() +
  labs(title="Tufte Styled Boxplot", 
       subtitle="Lower level SSA-SPVe",
       caption="Source: mpg",
       x="Node degreee",
       y="Sweetpotato regions")

#----- species.strength
ggplot(nohl, aes(IDs, species.strength))+
  geom_tufteboxplot(voffset = 0.0, hoffset = 0, median.type = "line") + 
  # geom_point()+
  stat_summary(fun = "mean", geom = "point", col = "black", size = 3, shape = 1)+
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  # theme_tufte() +
  labs(title="Tufte Styled Boxplot", 
       subtitle="Lower level SSA-SPVe",
       caption="Source: mpg",
       x="Node degreee",
       y="Sweetpotato regions")

#----- species.specificity.index
ggplot(nohl, aes(IDs, species.specificity.index))+
  geom_tufteboxplot(voffset = 0.0, hoffset = 0, median.type = "line") + 
  # geom_point()+
  stat_summary(fun = "mean", geom = "point", col = "black", size = 3, shape = 1)+
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  # theme_tufte() +
  labs(title="Tufte Styled Boxplot", 
       subtitle="Lower level SSA-SPVe",
       caption="Source: mpg",
       x="Node degreee",
       y="Sweetpotato regions")

#----- weighted.betweenness
ggplot(nohl, aes(IDs, weighted.betweenness))+
  geom_tufteboxplot(voffset = 0.0, hoffset = 0, median.type = "line") + 
  # geom_point()+
  stat_summary(fun = "mean", geom = "point", col = "black", size = 3, shape = 1)+
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  # theme_tufte() +
  labs(title="Tufte Styled Boxplot", 
       subtitle="Lower level SSA-SPVe",
       caption="Source: mpg",
       x="Node degreee",
       y="Sweetpotato regions")

#----- weighted.closeness
ggplot(nohl, aes(IDs, weighted.closeness))+
  geom_tufteboxplot(voffset = 0.0, hoffset = 0, median.type = "line") + 
  # geom_point()+
  stat_summary(fun = "mean", geom = "point", col = "black", size = 3, shape = 1)+
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  # theme_tufte() +
  labs(title="Tufte Styled Boxplot", 
       subtitle="Lower level SSA-SPVe",
       caption="Source: mpg",
       x="Node degreee",
       y="Sweetpotato regions")

#----- Fisher.alpha
ggplot(nohl, aes(IDs, Fisher.alpha))+
  geom_tufteboxplot(voffset = 0.0, hoffset = 0, median.type = "line") + 
  # geom_point()+
  stat_summary(fun = "mean", geom = "point", col = "black", size = 3, shape = 1)+
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  # theme_tufte() +
  labs(title="Tufte Styled Boxplot", 
       subtitle="Lower level SSA-SPVe",
       caption="Source: mpg",
       x="Node degreee",
       y="Sweetpotato regions")

#----- Fisher.alpha
ggplot(nohl, aes(IDs, Fisher.alpha))+
  geom_tufteboxplot(voffset = 0.0, hoffset = 0, median.type = "line") + 
  # geom_point()+
  stat_summary(fun = "mean", geom = "point", col = "black", size = 3, shape = 1)+
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  # theme_tufte() +
  labs(title="Tufte Styled Boxplot", 
       subtitle="Lower level SSA-SPVe",
       caption="Source: mpg",
       x="Node degreee",
       y="Sweetpotato regions")

#----- Fisher.alpha
ggplot(nohl, aes(IDs, Fisher.alpha))+
  geom_tufteboxplot(voffset = 0.0, hoffset = 0, median.type = "line") + 
  # geom_point()+
  stat_summary(fun = "mean", geom = "point", col = "black", size = 3, shape = 1)+
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  # theme_tufte() +
  labs(title="Tufte Styled Boxplot", 
       subtitle="Lower level SSA-SPVe",
       caption="Source: mpg",
       x="Node degreee",
       y="Sweetpotato regions")

#----- Fisher.alpha
ggplot(nohl, aes(IDs, Fisher.alpha))+
  geom_tufteboxplot(voffset = 0.0, hoffset = 0, median.type = "line") + 
  # geom_point()+
  stat_summary(fun = "mean", geom = "point", col = "black", size = 3, shape = 1)+
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  # theme_tufte() +
  labs(title="Tufte Styled Boxplot", 
       subtitle="Lower level SSA-SPVe",
       caption="Source: mpg",
       x="Node degreee",
       y="Sweetpotato regions")




#------ Node metrics LL
files3 = list.files(pattern="*metrics_ll.csv")
names <- strsplit(files3, split = "_")
#
for (i in 1:length(files3)){
  x <- read.csv(files3[i], as.is=T)
  x$IDs <- rep(names[[i]][1], nrow(x))
  ls[[i]] = x
  colnames(ls[[i]]) <- colnames(ls[[1]])
}

noll <- rbind(ls[[1]], ls[[2]], ls[[3]], ls[[4]], ls[[5]], ls[[6]], ls[[7]])
as_tibble(noll)

#===
ggplot(noll, aes(IDs, degree))+ 
  geom_count(aes(group = IDs),shape = 1, show.legend=T) +
  # scale_size_area()+
  labs(subtitle="Higher level SSA-SPV", 
       y="Node degree", 
       x="Sweetpotato region", 
       title="Count Points")


#------
library(ggthemes)
# plot
ggplot(noll, aes(IDs, degree))+
  geom_tufteboxplot(voffset = 0.0, hoffset = 0, median.type = "line") + 
  # geom_point()+
  stat_summary(fun = "mean", geom = "point", col = "black", size = 3, shape = 1)+
  theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
  # theme_tufte() +
  labs(title="Tufte Styled Boxplot", 
       subtitle="Lower level SSA-SPVe",
       caption="Source: mpg",
       x="Node degreee",
       y="Sweetpotato regions")










#----------------

hl1 <- read.csv("k-cluster1_node-metrics_hl.csv", as.is=T)
hl1

ll1 <- read.csv("k-cluster1_node-metrics_ll.csv", as.is=T) 
ll1
mean(ll1$degree)
max(ll1$degree)
