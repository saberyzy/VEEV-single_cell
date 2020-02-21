# usage: visualize genes in DE analysis is not affected by time
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

#plotting
VEEV_2 <- SetAllIdent(VEEV_2, id = 'viral_load_ridge')
VEEV_uninfect <- SubsetData(VEEV_2, ident.use = 'uninfect')
VEEV_uninfect <- SetAllIdent(VEEV_uninfect, id = 'hours')
names <- c('TNFAIP3', 'TAF7', 'RND3')
Vln_hours <- function(names){
  plist <- list()
  for (x in names){
    plist[[x]] <- VlnPlot(object = VEEV_uninfect, features.plot = x, point.size.use = FALSE) + 
      labs(x = '') +
      geom_boxplot(width=0.1, color="black", alpha=0.7, outlier.shape = NA) +
      theme_classic() +
      theme( 
        axis.text.x = element_text(face="bold", 
                                   size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(face="bold", 
                                   size=14),
        legend.position = "none")
  }
  return(plist)
}

plist <- Vln_hours(names)
p <- grid.arrange(grobs = plist, ncol = 1)
