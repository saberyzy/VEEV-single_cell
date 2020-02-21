# usage: plotting quality controls
# date updated: 02/19/2020
# author: Zhiyuan Yao


library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

#preprocessing data
VEEV.data <- read.table('/Users/yzhiyuan/workspace/VEEV/counts_process/VEEVcounts.tsv', sep='\t', row.names=1, header=T)
VEEV <- CreateSeuratObject(raw.data = VEEV.data, project = "VEEV")
ERCC <- grep(pattern = "ERCC-", x = rownames(x = VEEV@data), value = TRUE)
percent.ERCC <- Matrix::colSums(VEEV@raw.data[ERCC, ])/Matrix::colSums(VEEV@raw.data)
VEEV <- AddMetaData(object = VEEV, metadata = percent.ERCC, col.name = "percent.ERCC")
VEEV@meta.data["Total_reads"] <- VEEV@meta.data['nUMI']

#Total reads
VlnPlot(object = VEEV, features.plot = c("Total_reads")) +
  labs(x = '', y = 'Reads number', title = 'Total reads') +
  theme_classic() +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 15), 
        axis.text.x = element_text(face="bold", 
                                   size=14),
        axis.text.y = element_text(face="bold", 
                                   size=14), 
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size=16,face="bold")
  ) +
  scale_y_continuous(breaks = trans_breaks("log10", function(x) 10^x),
                     trans = log10_trans(),
                     labels = trans_format("log10", math_format(10^.x)))

#number of genes
VlnPlot(object = VEEV, features.plot = c("nGene")) +
  labs(x = '', y = 'Gene counts', title = 'Number of genes') +
  theme_classic() +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 15), 
        axis.text.x = element_text(face="bold", 
                                   size=14),
        axis.text.y = element_text(face="bold", 
                                   size=14), 
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size=16,face="bold")
  ) 

#ERCC percentage
VlnPlot(object = VEEV, features.plot = c("percent.ERCC")) +
  labs(x = '', y = 'ERCC ratio', title = 'ERCC ratio') +
  theme_classic() +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 15), 
        axis.text.x = element_text(face="bold", 
                                   size=14),
        axis.text.y = element_text(face="bold", 
                                   size=14), 
        axis.title=element_text(size=14,face="bold"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size=16,face="bold")
  )

VEEV@meta.data$totalident <- 'VEEVsample'
VEEV <- SetAllIdent(VEEV, id = "totalident")
