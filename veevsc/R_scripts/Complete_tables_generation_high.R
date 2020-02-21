# usage: Gene differential expression analysis between viral_load_high and uninfect,
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

# find differential expressed genes in viral_load_high_group
VEEV_2 <- SetAllIdent(VEEV_2, id = "viral_load_ridge")
viralloadmarker_high <- FindMarkers(VEEV_2, ident.1 = "viral_load_high", ident.2 = "uninfect", 
                                    min.pct = 0.3)
viralloadmarker_high_2 <- viralloadmarker_high[viralloadmarker_high$p_val_adj < 0.05, ]

# batch/time effect correction by correlation
# in uninfect cells
VEEV_2 <- SetAllIdent(VEEV_2, id = "MOI")
uninfectcells <- row.names(subset(metadata, viral_load_ridge == 'uninfect'))
VEEV_uninfect_data <- VEEV_2_nordata[,uninfectcells]
metadata_uninfect <- metadata[metadata[,'MOI'] == 0,]
metadata_uninfect[, 'hours_number'] <- 0
metadata_uninfect[ metadata_uninfect[,'hours'] == '0.5h', 'hours_number'] <- 0.5
metadata_uninfect[ metadata_uninfect[,'hours'] == '1.5h', 'hours_number'] <- 1.5
metadata_uninfect[ metadata_uninfect[,'hours'] == '4h', 'hours_number'] <- 4
metadata_uninfect[ metadata_uninfect[,'hours'] == '6h', 'hours_number'] <- 6
metadata_uninfect[ metadata_uninfect[,'hours'] == '12h', 'hours_number'] <- 12
metadata_uninfect[ metadata_uninfect[,'hours'] == '24h', 'hours_number'] <- 24
q <- as.numeric(as.character(metadata_uninfect$hours_number))
VEEV_uninfect_data['hours_number',] <- q

rownames_uninfect <- row.names(VEEV_uninfect_data)
pearson_uninfect <- data.frame(matrix(nrow=length(rownames_uninfect), ncol=2), row.names = rownames_uninfect)
colnames(pearson_uninfect) <- c('pearson_cor', 'p-value')

VEEV_uninfect_data_t <- as.data.frame(t(VEEV_uninfect_data))

for (gene in rownames_uninfect){
  targene <- VEEV_uninfect_data_t[,gene]
  testcor<- cor.test(targene, VEEV_uninfect_data_t$hours_number, method="pearson")
  pearson_uninfect[gene, 'pearson_cor'] <- testcor$estimate
  pearson_uninfect[gene, 'p-value'] <- testcor$p.value
}

pearson_uninfect[is.na(pearson_uninfect)] <- 0
pearson_uninfect_sort <- pearson_uninfect[order(pearson_uninfect$pearson_cor),]
pearson_uninfect_sort$adjust_p_value <- (pearson_uninfect_sort$`p-value`) * (length(rownames_uninfect))
pearson_uninfect_sort <- pearson_uninfect_sort[pearson_uninfect_sort[,'adjust_p_value'] < 0.05,]
pearson_uninfect_sort <- pearson_uninfect_sort[pearson_uninfect_sort[,'pearson_cor'] != 0,]


# batch/time effect correction by correlation
# in infected cells

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

rownames_infect <- row.names(VEEV_infect_data)
pearson_infect <- data.frame(matrix(nrow=length(rownames_infect), ncol=2), row.names = rownames_infect)
colnames(pearson_infect) <- c('pearson_cor', 'p-value')

VEEV_infect_data_t <- as.data.frame(t(VEEV_infect_data))
for (gene in rownames_infect){
  targene <- VEEV_infect_data_t[,gene]
  testcor<- cor.test(targene, VEEV_infect_data_t$hours_number, method="pearson")
  pearson_infect[gene, 'pearson_cor'] <- testcor$estimate
  pearson_infect[gene, 'p-value'] <- testcor$p.value
}

pearson_infect[is.na(pearson_infect)] <- 0
pearson_infect_sort <- pearson_infect[order(pearson_infect$pearson_cor),]
pearson_infect_sort$adjust_p_value <- (pearson_infect_sort$`p-value`) * (length(rownames_infect))
pearson_infect_sort <- pearson_infect_sort[pearson_infect_sort['adjust_p_value'] < 0.05,]
pearson_infect_sort <- pearson_infect_sort[pearson_infect_sort[,'pearson_cor'] != 0,]

# generate a complete table for future use 
spearman_ercc_exp = read.table("/Users/yzhiyuan/workspace/VEEV/VEEV_Rproject/table/spearman_ercc.tsv",
                               sep='\t', row.names=1, header=T)
colnames(viralloadmarker_high_2) <- c('p_val_high', 'avg_logFC_high', 'pct.1_high', 'pct.2_high', 'p_val_adj_high')
a <- viralloadmarker_high_2
b <- intersect(row.names(a), row.names(spearman_ercc_exp))
c <- spearman_ercc_exp[b,]
a <- a[b, ]
temp <- cbind(a, c)
colnames(temp)[6] <- 'spearman_cor'
d <- intersect(row.names(temp), row.names(pearson_infect))
f <- pearson_infect[d,]
colnames(f) <- c('pearson_cor_infect', 'p_val_pearson_infect')
temp <- temp[d,]
temp <- cbind(temp, f)
g <- pearson_uninfect[d,]
colnames(g) <- c('pearson_cor_uninfect', 'p_val_pearson_uninfect')
temp <- cbind(temp, g)
Completelist_high <- temp


write.table(Completelist_high, file = '/Users/yzhiyuan/workspace/VEEV/VEEV_Rproject/table/Completelist_high.tsv', 
            quote = FALSE, sep = '\t')





