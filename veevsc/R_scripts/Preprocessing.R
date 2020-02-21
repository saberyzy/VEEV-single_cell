# usage: preprocessing data from counts table
# date updated: 02/19/2020
# author: Zhiyuan Yao

library(Seurat)
library(plyr)
library(dplyr)
library(reshape2)
library(scales)
library(ggrepel)
library(tidyverse)
library(ggplot2)

# Load data and check the quality 
VEEV.data <- read.table('/Users/yzhiyuan/workspace/VEEV/counts_process/VEEVcounts.tsv', sep='\t', row.names=1, header=T)
VEEV <- CreateSeuratObject(raw.data = VEEV.data, project = "VEEV")
ERCC <- grep(pattern = "ERCC-", x = rownames(x = VEEV@data), value = TRUE)
percent.ERCC <- Matrix::colSums(VEEV@raw.data[ERCC, ])/Matrix::colSums(VEEV@raw.data)
VEEV <- AddMetaData(object = VEEV, metadata = percent.ERCC, col.name = "percent.ERCC")
VEEV@meta.data["Total_reads"] <- VEEV@meta.data['nUMI']
VlnPlot(object = VEEV, features.plot = c("Total_reads"))
VlnPlot(object = VEEV, features.plot = c("nGene", "percent.ERCC"), nCol = 2, return = TRUE)
VEEV@meta.data$totalident <- 'VEEVsample'
VEEV <- SetAllIdent(VEEV, id = "totalident")

# make the metatable, labeling cells with MOI and hours

metadata <- VEEV@meta.data
a <- metadata[grepl('1001703003', VEEV@meta.data$orig.ident), ]
metadata$hours <- ''

# label by hours
metadata$hours <- ifelse(grepl('1001703001', metadata$orig.ident), '0.5h', 
                         ifelse(grepl('1001703003', metadata$orig.ident), '1.5h',
                                ifelse(grepl('1001703005', metadata$orig.ident),'4h',
                                       ifelse(grepl('1001703007', metadata$orig.ident), '6h', 
                                              ifelse(grepl('1001703009', metadata$orig.ident), '12h', 
                                                     '24h')))))

metadata$hours <- factor(metadata$hours, levels = c('0.5h', '1.5h', '4h', '6h', '12h', '24h'))
VEEV@meta.data <-metadata

#label by MOIs
metadata_new <- VEEV@meta.data
metadata_new$MOI <- ''
metadata_new[grep(pattern = "1001703001_[A-Z][1-6]$", rownames(metadata_new), value = TRUE), 'MOI'] <- 0
metadata_new[grep(pattern = "1001703001_[A-Z][7-9]$", rownames(metadata_new), value = TRUE), 'MOI'] <- 0.1
metadata_new[grep(pattern = "1001703001_[A-Z]1[0-2]$", rownames(metadata_new), value = TRUE), 'MOI'] <- 0.1
metadata_new[grep(pattern = "1001703001_[A-Z]1[3-8]$", rownames(metadata_new), value = TRUE), 'MOI'] <- 1
metadata_new[grep(pattern = "1001703001_[A-Z]19$", rownames(metadata_new), value = TRUE), 'MOI'] <- 3
metadata_new[grep(pattern = "1001703001_[A-Z]2[0-4]$", rownames(metadata_new), value = TRUE), 'MOI'] <- 3

metadata_new[grep(pattern = "10017030[0-1][2-9]_[A-Z][1-8]$", rownames(metadata_new), value = TRUE), 'MOI'] <- 0
metadata_new[grep(pattern = "10017030[0-1][2-9]_[A-Z]9$", rownames(metadata_new), value = TRUE), 'MOI'] <- 0.1
metadata_new[grep(pattern = "10017030[0-1][2-9]_[A-Z]1[0-6]$", rownames(metadata_new), value = TRUE), 'MOI'] <- 0.1
metadata_new[grep(pattern = "10017030[0-1][2-9]_[A-Z]1[7-9]$", rownames(metadata_new), value = TRUE), 'MOI'] <- 1
metadata_new[grep(pattern = "10017030[0-1][2-9]_[A-Z]2[0-4]$", rownames(metadata_new), value = TRUE), 'MOI'] <- 1

metadata_new[grep(pattern = "1001703011_[A-Z][1-8]$", rownames(metadata_new), value = TRUE), 'MOI'] <- 0
metadata_new[grep(pattern = "1001703011_[A-Z]9$", rownames(metadata_new), value = TRUE), 'MOI'] <- 0.1
metadata_new[grep(pattern = "1001703011_[A-Z]1[0-6]$", rownames(metadata_new), value = TRUE), 'MOI'] <- 0.1
metadata_new[grep(pattern = "1001703011_[A-Z]1[7-9]$", rownames(metadata_new), value = TRUE), 'MOI'] <- 1
metadata_new[grep(pattern = "1001703011_[A-Z]2[0-4]$", rownames(metadata_new), value = TRUE), 'MOI'] <- 1

metadata_new$MOI <- factor(metadata_new$MOI, level = c(0, 0.1, 1, 3))
VEEV@meta.data <-metadata_new

#filter cells
VEEV <- FilterCells(object = VEEV, subset.names = c("nGene", "percent.ERCC", "nUMI"), 
                    low.thresholds = c(500, -Inf, 300000), high.thresholds = c(Inf, 0.05, Inf))

#filter genes
MITO <- grep(pattern = "MT-", x = rownames(x = VEEV@data), value = TRUE)
RPS <- grep(pattern = "RPS", x = rownames(x = VEEV@data), value = TRUE)
RPL <- grep(pattern = "RPL", x = rownames(x = VEEV@data), value = TRUE)
SP <- grep(pattern = "RNA5S", x = rownames(x = VEEV@data), value = TRUE)
others <- grep(pattern = "__", x = rownames(x = VEEV@data), value = TRUE)
others <- c(others, ERCC)
non_nuclear <- c(MITO, RPS, RPL, SP)
excludeall <- c(others, non_nuclear)
VEEV_2 <- SubsetData(object = VEEV, ident.use = 'VEEVsample')
cell_metadata <- row.names(VEEV_2@meta.data)
VEEV_ercc <- VEEV@raw.data[, cell_metadata]
VEEV_ercc <- VEEV_ercc[!row.names(VEEV_ercc) %in% excludeall, ]
metadata <- VEEV_2@meta.data

#add ERCC and host reads information to metadata
ERCC_data <- VEEV@raw.data[ERCC, cell_metadata]
ERCC_sum <- data.frame(colname = names(ERCC_data),colSums_ERCC=colSums(ERCC_data))
metadata$ERCC_sum <- ERCC_sum$colSums_ERCC



#ERCC based normalization
a  <- VEEV_ercc / (as.vector(metadata$ERCC_sum))
a <- log((a*10000)  +  1)
b <- as.matrix(a)
c <- as(b, "sparseMatrix")
VEEV_2@data <- c
VEEV_2_nordata <- as.data.frame(x=as.matrix(x=VEEV_2@data))

#add virus read information to the metadata
metadata <-VEEV_2@meta.data
metadata$virus_read <- 0
a <- VEEV_nordata['VEEVGFP',VEEV_nordata['VEEVGFP',] != 0]
colname <- colnames(a)
for(name in colname){
  metadata[name, 'virus_read'] <- VEEV_2_nordata['VEEVGFP', name]
}
metadata$virus_read_raw <- 0
for(name in colname){
  metadata[name, 'virus_read_raw'] <- VEEV_2@raw.data['VEEVGFP', name]
}
metadata$Total_ERCC_ratio <- (metadata$Total_reads - metadata$ERCC_sum - metadata$virus_read_raw) / metadata$ERCC_sum
metadata$Host_reads <- metadata$Total_reads - metadata$ERCC_sum - metadata$virus_read_raw
metadata$virus_ERCC_ratio <- metadata$virus_read_raw / metadata$ERCC_sum
metadata$Virus_Total_ratio <- metadata$virus_read_raw / metadata$Total_reads
metadata$Host_Total_ratio <- metadata$Host_reads / metadata$Total_reads

metadata$infect <- 'No'
a <- rownames(subset(metadata, virus_read_raw >= 10, ))
metadata[a, 'infect'] <- 'Yes'

#save
VEEV_2@meta.data <- metadata
saveRDS(VEEV_2, file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")



