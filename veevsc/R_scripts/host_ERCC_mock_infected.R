# usage: Host_ERCC ratio in mock-infected cells
# date updated: 02/19/2020
# author: Zhiyuan Yao
# other scripts needed: Preprocessing.R

library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

# load Rproject
VEEV_2 <- readRDS(file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")
metadata <- VEEV_2@meta.data
metadata$host_ERCC_ratio <- (metadata$Host_reads) / (metadata$ERCC_sum)
metadata_uninfect <- subset(metadata, MOI == 0, )
a <- metadata_uninfect[, c('hours', 'host_ERCC_ratio')]

# plotting
p <- ggplot(a, aes(hours,host_ERCC_ratio, fill = hours)) + 
  geom_violin(trim = FALSE) +
  geom_jitter(shape=16, position=position_jitter(0.2))

p +
  labs(title='Host/ERCC ratio in uninfected cells', y="Host/ERCC reads ratio", x = "Hours") +
  theme_classic() +
  theme(legend.title = element_text(size = 20),
      legend.text = element_text(size = 15), 
      axis.text.x = element_text(face="bold", 
                                 size=14, angle = 45, hjust = 1),
      axis.text.y = element_text(face="bold", 
                                 size=14), 
      axis.title=element_text(size=14,face="bold") 
)
