#'@title:  "Histogram of sweetpotato virome: cluster"
#'@author: "R Alcala"
#'@date:   "02/02/2021"
#'@output: "Metadata"

# library
library(ggplot2)
library(tidytext)

# ViNA results data 
virome <- read.csv("papacoldland-bipartitenetworkMar15.RData", header = T, stringsAsFactors = F)[-1]
# virome <- merge(virome, kclusters[c(2,8)], 
#                 by.x= "IDs", by.y= "SampleID")

as_tibble(virome)

# Barplot sweetpotato clusters stacked + percent
pdf(paste0("papa_virome_altzones_barplot_regions_",format(Sys.time(), "%b%d"), ".pdf"),
    width = 30, # The width of the plot in inches
    height = 15) # The height of the plot in inches
  ggplot(virome, aes(x=reorder(IDs, Family), y=log2(RPKM_mean), fill=Genus)) + 
    geom_bar(position="fill", stat="identity", width = 1, alpha = 0.75)+
    ggtitle("Virome genus frequency by sweetpotato land area ") +
    facet_grid(~altzones, scales = "free", space = "free") +
    scale_fill_viridis(discrete = T, option = "B") +
    scale_x_reordered() +
    theme_bw()
dev.off()
