# usage: unsupervised clustering for the data and visualize the clusters by TSNE 
# date updated: 02/19/2020
# author: Zhiyuan Yao
# other scripts needed : Preprocessing.R

library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

# load data
VEEV_2 <- readRDS(file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")
metadata <- VEEV_2@meta.data
VEEV_2 <-SetAllIdent(VEEV_2, id = 'MOI')

# select infected cells
infectcells <- row.names(subset(metadata, MOI != 0, ))
uninfectcells <- row.names(subset(metadata, MOI == 0, ))
VEEV_cluster <- SubsetData(object = VEEV_2, cells.use = infectcells)

# clustering
VEEV_cluster <- FindVariableGenes(object = VEEV_cluster, mean.function = ExpMean, dispersion.function = LogVMR, 
                                  x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 2)
length(VEEV_cluster@var.genes)
VEEV_cluster <- ScaleData(object = VEEV_cluster)
VEEV_cluster <- RunPCA(object = VEEV_cluster, pc.genes = VEEV_cluster@var.genes)
PCElbowPlot(object = VEEV_cluster)
VEEV_cluster <- FindClusters(object = VEEV_cluster, reduction.type = "pca", dims.use = 1:15, 
                             resolution = 0.6, print.output = 0, save.SNN = TRUE)
VEEV_cluster <- RunTSNE(object = VEEV_cluster, dims.use = 1:15, do.fast = TRUE)

# visualization 
TSNEPlot(object = VEEV_cluster)
VlnPlot(object = VEEV_cluster, features.plot = c("VEEVGFP"))
