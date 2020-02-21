# usage: Gene differential expression analysis between viral_load_low and uninfect,
#        generate the information table for future use
# date updated: 02/19/2020
# author: Zhiyuan Yao
# other scripts needed : grouping_by_viral_load.R, spearman_histo.R

library(Seurat)
library(plyr)
library(dplyr)

# Load Rprojects
VEEV_2 <- readRDS(file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")
metadata <- VEEV_2@meta.data
VEEV_2_nordata <- as.data.frame(x=as.matrix(x=VEEV_2@data))

# find differential expressed genes in viral_load_low_group
VEEV_2 <- SetAllIdent(VEEV_2, id = "viral_load_ridge")
viralloadmarker_low <- FindMarkers(VEEV_2, ident.1 = "viral_load_low", ident.2 = "uninfect", 
                                    min.pct = 0.3)
viralloadmarker_low_2 <- viralloadmarker_low[viralloadmarker_low$p_val_adj < 0.05, ]

# time/batch correction by correlation in uninfected cells
uninfectcells <- row.names(metadata[metadata[,'MOI'] == 0,])
metadata_uninfect <- metadata[metadata[,'MOI'] == 0,]
VEEV_uninfect_data <- VEEV_2_nordata[,uninfectcells]
metadata_uninfect[, 'hours_number'] <- 0

metadata_uninfect[ metadata_uninfect[,'hours'] == '0.5h', 'hours_number'] <- 0.5
metadata_uninfect[ metadata_uninfect[,'hours'] == '1.5h', 'hours_number'] <- 1.5
metadata_uninfect[ metadata_uninfect[,'hours'] == '4h', 'hours_number'] <- 4
metadata_uninfect[ metadata_uninfect[,'hours'] == '6h', 'hours_number'] <- 6
metadata_uninfect[ metadata_uninfect[,'hours'] == '12h', 'hours_number'] <- 12
metadata_uninfect[ metadata_uninfect[,'hours'] == '24h', 'hours_number'] <- 24

q <- as.numeric(as.character(metadata_uninfect$hours_number))
VEEV_uninfect_data['hours_number',] <- q

h <- c(0.5, 1.5, 4, 6)
VEEV_uninfect_data_early <- VEEV_uninfect_data[, VEEV_uninfect_data['hours_number',] %in% h]

rownames_uninfect_early <- row.names(VEEV_uninfect_data_early)
pearson_uninfect_early <- data.frame(matrix(nrow=length(rownames_uninfect_early), ncol=2), 
                                     row.names = rownames_uninfect_early)
colnames(pearson_uninfect_early) <- c('pearson_cor', 'p-value')

VEEV_uninfect_data_early_t <- as.data.frame(t(VEEV_uninfect_data_early))

for (gene in rownames_uninfect_early){
  targene <- VEEV_uninfect_data_early_t[,gene]
  testcor<- cor.test(targene, VEEV_uninfect_data_early_t$hours_number, method="pearson")
  pearson_uninfect_early[gene, 'pearson_cor'] <- testcor$estimate
  pearson_uninfect_early[gene, 'p-value'] <- testcor$p.value
}

pearson_uninfect_early[is.na(pearson_uninfect_early)] <- 0

# time/batch correction by correlation in infected cells 
metadata_infect <- metadata[metadata[,'MOI'] != 0,]
metadata_infect[, 'hours_number'] <- 0
VEEV_infect_data <- VEEV_2_nordata[,row.names(metadata_infect)]

metadata_infect[ metadata_infect[,'hours'] == '0.5h', 'hours_number'] <- 0.5
metadata_infect[ metadata_infect[,'hours'] == '1.5h', 'hours_number'] <- 1.5
metadata_infect[ metadata_infect[,'hours'] == '4h', 'hours_number'] <- 4
metadata_infect[ metadata_infect[,'hours'] == '6h', 'hours_number'] <- 6
metadata_infect[ metadata_infect[,'hours'] == '12h', 'hours_number'] <- 12
metadata_infect[ metadata_infect[,'hours'] == '24h', 'hours_number'] <- 24

q <- as.numeric(as.character(metadata_infect$hours_number))
VEEV_infect_data['hours_number',] <- q

VEEV_infect_data_early <- VEEV_infect_data[, VEEV_infect_data['hours_number',] %in% h]

rownames_infect_early <- row.names(VEEV_infect_data_early)
pearson_infect_early <- data.frame(matrix(nrow=length(rownames_infect_early), ncol=2), 
                                   row.names = rownames_infect_early)
colnames(pearson_infect_early) <- c('pearson_cor', 'p-value')

VEEV_infect_data_early_t <- as.data.frame(t(VEEV_infect_data_early))

for (gene in rownames_infect_early){
  targene <- VEEV_infect_data_early_t[,gene]
  testcor<- cor.test(targene, VEEV_infect_data_early_t$hours_number, method="pearson")
  pearson_infect_early[gene, 'pearson_cor'] <- testcor$estimate
  pearson_infect_early[gene, 'p-value'] <- testcor$p.value
}

pearson_infect_early[is.na(pearson_infect_early)] <- 0

# generate a complete table for future use 
spearman_ercc_exp = read.table("/Users/yzhiyuan/workspace/VEEV/VEEV_Rproject/table/spearman_ercc.tsv",
                               sep='\t', row.names=1, header=T)
colnames(viralloadmarker_low_2) <- c('p_val_low', 'avg_logFC_low', 'pct.1_low', 'pct.2_low', 'p_val_adj_low')
a <- viralloadmarker_low_2
b <- intersect(row.names(a), row.names(spearman_ercc_exp))
c <- spearman_ercc_exp[b,]
a <- a[b, ]
temp <- cbind(a, c)
colnames(temp)[6] <- 'spearman_cor'
d <- intersect(row.names(temp), row.names(pearson_infect_early))
f <- pearson_infect_early[d,]
colnames(f) <- c('pearson_cor_infect_early', 'p_val_pearson_infect_early')
temp <- temp[d,]
temp <- cbind(temp, f)
g <- pearson_uninfect_early[d,]
colnames(g) <- c('pearson_cor_uninfect_early', 'p_val_pearson_uninfect_early')
temp <- cbind(temp, g)
Completelist_low <- temp

write.table(Completelist_low, file = '/Users/yzhiyuan/workspace/VEEV/VEEV_Rproject/table/Completelist_low_new.tsv', 
            quote = FALSE, sep = '\t')








