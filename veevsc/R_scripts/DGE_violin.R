# usage: visualize the gene differential expression of several representative genes
# date updated: 02/19/2020
# author: Zhiyuan Yao
# other scripts needed : grouping_by_viral_load.R

library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)

#load data
VEEV_2 <- readRDS(file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")
metadata <- VEEV_2@meta.data
VEEV_2_nordata <- as.data.frame(x=as.matrix(x=VEEV_2@data))

# plotting
VEEV_2 <- SetAllIdent(VEEV_2, id = 'viral_load_ridge')
VEEV_2@ident <- factor(x = VEEV_2@ident, levels = c('uninfect', 'viral_load_low', 'viral_load_high'))
names <- c('TNFAIP3', 'TAF7', 'RND3')
Vln <- function(names){
  plist <- list()
  for (x in names){
    plist[[x]] <- VlnPlot(object = VEEV_2, features.plot = x, point.size.use = FALSE) + 
      labs(x = '') +
      geom_boxplot(width=0.1, color="black", alpha=0.7, outlier.shape = NA) +
      theme( 
        axis.text.x = element_text(face="bold", 
                                   size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(face="bold", 
                                   size=14)
      )
  }
  return(plist)
}

plist <- Vln(names)
p <- grid.arrange(grobs = plist, ncol = length(plist))


