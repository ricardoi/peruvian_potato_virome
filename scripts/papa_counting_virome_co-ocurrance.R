# Code to species co-occurrence
#### NOTE
# k1vina comes from the script ViNAbn_v2.hpg.R

viromeReads <- ddply(virome2, .(IDs, Family, Genus, pSpecies, Acronym, altzones), 
                     summarise, Coverage=mean(Reads_mean))

virome0 <- as_tibble(ddply(virome2, .(IDs, Family, Genus, pSpecies, Acronym, altzones), summarise, Coverage=mean(RPKM_mean)))
virome0

#virome0=viromeReads[viromeReads$Coverage > 1,]
IDs.mastre <- virome0$IDs[virome0$Acronym == "SPSMV-1"]
mastre.0 <- virome0[which(virome0$IDs %in% IDs.mastre), ]
table(mastre.0$Acronym)

bin = table(virome0$IDs, virome0$Acronym)


pair1 <- bin[,c("SPFMV", "SPCSV", "SPLCV", "SPSMV-1")]
pair1 <- data.frame(SPFMV=pair1[,1], SPCSV=pair1[,2], SPLVC=pair1[,3], SPSMV=pair1[,4])
pair1$mixed <- apply(pair1, 1, sum)
table(pair1$mixed)

distribution <- sort(table(virome0$Acronym), decreasing = T)
distribution[1] <- distribution[1]/2
plot(distribution, las = 2)

pair1 <- bin[,c("SPFMV", "SPCSV", "SPLCV", "SPSMV-1", "BYMV")]
pair1 <- data.frame(SPFMV=pair1[,1], SPCSV=pair1[,2], SPLVC=pair1[,3], SPSMV=pair1[,4], BYMV=pair1[,5])
pair1$mixed <- apply(pair1, 1, sum)
table(pair1$mixed)


pair1 <- bin[,c("SPPV", "CmV", "SPFMV", "SPSMV-1", "SPCSV", "SPVG", "SPLCV", "SPV2", "SPVC")]

pair1 <- data.frame(SPPV=pair1[,1], CmV=pair1[,2], SPFMV=pair1[,3], SPSMV=pair1[,4], SPCSV=pair1[,5],
                    SPLVC=pair1[,6], SPV2=pair1[,7], SPVC=pair1[,8])
pair1$mixed <- apply(pair1, 1, sum)
table(pair1$mixed)

