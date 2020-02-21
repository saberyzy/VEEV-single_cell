# usage: 3p/5p ratio correlation with genes
# date updated: 02/20/2020
# author: Zhiyuan Yao

library(Seurat)
library(plyr)
library(dplyr)
library(reshape2)
library(scales)
library(ggrepel)
library(tidyverse)
library(ggplot2)
library(gridExtra)

# load data
vRNA53<- read.table('/Users/yzhiyuan/workspace/VEEV/VEEV_Rproject/table/vRNA_5p_3p.tsv', sep='\t', 
                    row.names=1, header=T)
VEEV_2 <- readRDS(file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")
VEEV_2_nordata <- as.data.frame(x=as.matrix(x=VEEV_2@data))
metadata <- VEEV_2@meta.data

# Formalize data and Normalization 
a <- vRNA53
raw <- vRNA53
row.names(a) <- gsub(x = row.names(a), pattern = "10017", replacement = "X10017")
row.names(raw) <- gsub(x = row.names(raw), pattern = "10017", replacement = "X10017")
a <- log(a+1)
infectcells <- row.names(subset(metadata, MOI != 0, ))
vRNA_name <- row.names(a)
cell_vRNA <- intersect(infectcells, vRNA_name)
raw <- raw[cell_vRNA,]
d <- as.data.frame(cell_vRNA)
d$cell_vRNA <- gsub(x = d$cell_vRNA, pattern = "X10017", replacement = "10017")
b <- a[cell_vRNA,]
b$ratio_3p_5p <- (b$vRNA_3p) / (b$vRNA_5p)
b$raw_ratio_3p_5p <- (raw$vRNA_3p) / (raw$vRNA_5p)
raw$raw_ratio_3p_5p <- (raw$vRNA_3p) / (raw$vRNA_5p)


# try calculate spearman between gene expression level and the 3p/5p ratio
b_2 <- b[b$vRNA_log > 1,]
name_b_2 <- row.names(b_2)
VEEV_53 <- VEEV_2_nordata[,name_b_2]
VEEV_53_t <- as.data.frame(t(VEEV_53))
VEEV_53_t[,'vRNA_5p'] <- b_2$vRNA_5p
VEEV_53_t[,'vRNA_3p'] <- b_2$vRNA_3p
VEEV_53_t[,'raw_ratio_3p_5p'] <- b_2$raw_ratio_3p_5p

colnames_vRNA53 <- colnames(VEEV_53_t)
spearman_vRNA53 <- data.frame(matrix(nrow=length(colnames_vRNA53), ncol=3), row.names = colnames_vRNA53)
colnames(spearman_vRNA53) <- c('spearman_cor_3p', 'spearman_cor_5p', 'spearman_ratio_3p_5p')


for (gene in colnames_vRNA53){
  targene <- VEEV_53_t[,gene]
  testcor_5p <- cor.test(targene, VEEV_53_t$vRNA_5p, method="spearman", exact=FALSE)
  testcor_3p <- cor.test(targene, VEEV_53_t$vRNA_3p, method="spearman", exact=FALSE)
  testcor_ratio <- cor.test(targene, VEEV_53_t$raw_ratio_3p_5p, method="spearman", exact=FALSE)
  spearman_vRNA53[gene, 'spearman_cor_5p'] <- testcor_5p$estimate
  spearman_vRNA53[gene, 'spearman_cor_3p'] <- testcor_3p$estimate
  spearman_vRNA53[gene, 'spearman_ratio_3p_5p'] <- testcor_ratio$estimate
}

spearman_vRNA53[is.na(spearman_vRNA53)] <- 0
spearman_vRNA53_nozero <- subset(spearman_vRNA53, 
                                 (spearman_cor_5p != 0) & (spearman_cor_3p != 0) & (spearman_ratio_3p_5p != 0)
                                 , )
v = c('VEEVGFP', 'vRNA_5p', 'vRNA_3p', 'raw_ratio_3p_5p')
spearman_vRNA53_nozero <- subset(spearman_vRNA53_nozero , 
                                 !(row.names(spearman_vRNA53_nozero) %in% v),)

# plot histo for spearman_ratio_3p_5p
p <- ggplot(spearman_vRNA53_nozero, aes(x=spearman_ratio_3p_5p)) + 
  geom_histogram(binwidth=.025, colour="black", fill="white") 
p

# plot scatter 
genename <- c('PFN2', 'BROX', 'ATP6V1B2','BNIP3','LAMP2','PIP4K2A', 'VAMP7', 'RAB7A', 'SEC22B')
VEEV_2_ERCCnor  <- VEEV_2_nordata
VEEV_2_nordata_233 <- VEEV_2_ERCCnor[,name_b_2]
df <- VEEV_2_nordata_233[genename,]
df <- as.data.frame(t(df))
df$ratio_3p_5p <- b_2$raw_ratio_3p_5p

plist <- list()
Scatter <- function(genename){
  plist <- list()
  for (n in genename){
    plist[[n]] <- ggplot(df, aes_string(x='ratio_3p_5p', y=n)) + 
      stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon') + 
      scale_fill_continuous(low="green",high="red") +
      geom_point(alpha=0.3) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      labs(x = 'vRNA 3p reads / vRNA 5p reads', y = 'Gene reads / ERCC reads', title = n) +
      annotation_logticks(sides = "l") +
      theme_bw() +
      theme( 
        axis.text.x = element_text(face="bold", 
                                   size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(face="bold", 
                                   size=14),
        legend.position = 'none'
      )
    label <- paste('Rho = ',format(round(spearman_vRNA53_nozero[n,'spearman_ratio_3p_5p'], 2), nsmall = 2), sep='')
    print(label)
    plist[[n]] <- plist[[n]] + annotate("text", x = 2.5, y = 3, label = label)
  }
  return(plist)
}

plist <- Scatter(genename)
quartz()
p <- grid.arrange(grobs = plist, ncol = 3)
