# usage: compare superproducer cells and bystanders at 6h MOI == 1 
# date updated: 02/19/2020
# author: Zhiyuan Yao
# other scripts needed : Preprocessing.R

library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)

# load Rprojects
VEEV_2 <- readRDS(file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")
metadata <- VEEV_2@meta.data
VEEV_2_nordata <- as.data.frame(x=as.matrix(x=VEEV_2@data))

# normalized the data without log transfer in convenience of plotting
VEEV_ercc <- VEEV_2@raw.data[, row.names(metadata)]
a  <- VEEV_ercc / (as.vector(metadata$ERCC_sum))
VEEV_2_ERCCnor <- a[row.names(VEEV_2_nordata),]
VEEV_2_ERCCnor_2 <- VEEV_2_ERCCnor * 10000 + 1
m <- as.matrix(VEEV_2_ERCCnor_2)
c <- as(m, "sparseMatrix")
VEEV_2@data <- c

# choose subdata set to find target cells
cells <- row.names(subset(metadata, MOI == 1 & hours == '6h',))
VEEV_moi1_6h <- SubsetData(VEEV_2, cells.use = cells)
metadata_moi1_6h <- VEEV_moi1_6h@meta.data
metadata_moi1_6h_sort <- metadata_moi1_6h[order(metadata_moi1_6h$virus_read),]
cells_high <- row.names(tail(metadata_moi1_6h_sort, n = 13))
cells_low <- row.names(head(metadata_moi1_6h_sort, n = 13))
cells <- c(cells_high, cells_low)
rm(VEEV_moi1_6h)

# choose subdata for plotting
VEEV_small <- SubsetData(VEEV_2, cells.use = cells)
VEEV_small <- SetAllIdent(VEEV_small, id = 'viral_load_ridge')
Markers_high_low <- FindMarkers(VEEV_small, ident.1 = 'viral_load_high')
Markers_high_low_2 <- Markers_high_low[Markers_high_low$p_val < 0.05,]
Markers_high_low_2 <- Markers_high_low_2[order(Markers_high_low_2$avg_logFC, decreasing = TRUE),]
Markers_high_low$sbstract <- Markers_high_low$pct.1 - Markers_high_low$pct.2
genelist <- subset(Markers_high_low, sbstract >= 0.5 | sbstract <= -0.5,)
genelist <- genelist[order(genelist$sbstract, decreasing = TRUE),]

# plotting
names <- c('VEEVGFP', 'MEIS2', 'WWP1', 'ZMAT5', 'PARP1')
Vln_6h <- function(names){
  plist <- list()
  Seuratobject <- VEEV_small
  for (x in names){
    plist[[x]] <- VlnPlot(object = Seuratobject, features.plot = x) +
      scale_y_continuous(trans = 'log10', breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      labs(x = '', y = 'Gene reads / ERCC reads + 0.0001') +
      annotation_logticks(sides = "l", base = 10) +
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
  
plist <- Vln_6h(names) 
quartz()
p <- grid.arrange(grobs = plist, ncol = length(plist))

ggsave('/Users/yzhiyuan/workspace/VEEV/VEEV_scripts/figures/figure_1H/figure_1I_2.pdf', p, 
       width = 10, height = 4)



