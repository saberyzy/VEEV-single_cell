# usage: visualize the time/batch effect correction for viral_low_group
# date updated: 02/19/2020
# author: Zhiyuan Yao
# other scripts needed : Complete_tables_generation_low.R

library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

# load data
VEEV_2 <- readRDS(file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")
metadata <- VEEV_2@meta.data
VEEV_2_nordata <- as.data.frame(x=as.matrix(x=VEEV_2@data))
Completelist_low <- read.table('/Users/yzhiyuan/workspace/VEEV/VEEV_Rproject/table/Completelist_low_new.tsv', 
                               sep='\t', row.names=1, header=T) 

# mark genes not affected by time/batch effect 
a <- c('pearson_cor_infect_early', 'pearson_cor_uninfect_early')
b <- Completelist_low[,a]
b$pearson_diff <- b$pearson_cor_infect / b$pearson_cor_uninfect
b$name <- row.names(b)
name_labeldata_1 <- row.names(subset(b, pearson_diff > 2 | pearson_diff < 0.5,))
c <- row.names(b) %in% c(name_labeldata_1)
default <- row.names(b[!c,])
b[name_labeldata_1,'group'] <- 'Positive signal'
b[default,'group'] <- 'Batch effect'

# plotting
ggplot(b, aes(pearson_cor_uninfect_early, pearson_cor_infect_early)) + 
  geom_point(size = 2, alpha = 0.4, aes(color = group)) +
  geom_abline(intercept = 0, slope = 0.5, , size = 1, linetype = "dashed") +
  geom_abline(intercept = 0, slope = 2, size = 1, linetype = "dashed") +
  labs(x = 'Correlation in uninfected cells', y = 'Correlation in infected cells', 
       title="Time-gene expression correlation ") +
  xlim(-0.6, 0.4) +
  ylim(-0.6, 0.4) +
  coord_fixed(ratio = 1) +
  theme_classic() +
  scale_color_manual(values=c('#999999','red')) +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 15), 
        axis.text.x = element_text(face="bold", 
                                   size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(face="bold", 
                                   size=14), 
        axis.title=element_text(size=12,face="bold") 
  )
