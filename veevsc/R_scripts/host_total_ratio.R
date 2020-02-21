# usage: host-total ratio in infected cells
# date updated: 02/19/2020
# author: Zhiyuan Yao
# other scripts needed: Preprocessing.R


library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

# load Rprojects
VEEV_2 <- readRDS(file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")
metadata <- VEEV_2@meta.data
metadata_infect <- subset(metadata, infect == 'Yes', )
a <- metadata_infect[, c('hours', 'Host_Total_ratio')]

# plotting
p <- ggplot(a, aes(hours,Host_Total_ratio, fill = hours)) + 
  geom_boxplot(outlier.shape = NA) 

p + labs(title='Host/Total reads ratio in infected cells', y="Host/Total reads ratio", x = "Hours Post Infection") +
  theme_classic() +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 15), 
        axis.text.x = element_text(face="bold", 
                                   size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(face="bold", 
                                   size=14), 
        axis.title=element_text(size=14,face="bold") 
  )

