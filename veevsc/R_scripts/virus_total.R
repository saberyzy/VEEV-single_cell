# usage: plotting virus-total ratio in infected cells
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
metadata_infect <- subset(metadata, infect == 'Yes', )
a <- metadata_infect[, c('hours', 'Virus_Total_ratio')]
total <- metadata_infect[,c('totalident','Virus_Total_ratio')]
total$totalident <- 'total'
colnames(total) <- c('hours','Virus_Total_ratio')
a <- rbind(a, total)


# plotting
p <- ggplot(a, aes(hours,Virus_Total_ratio, fill = hours)) + 
  geom_violin(scale = 'width') +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+
  geom_boxplot(width=0.1, color="black", alpha=0.7, outlier.shape = NA, fill = "white") +
  geom_hline(yintercept = 0.001, color = 'coral', linetype = 2) +
  annotation_logticks(sides = "l")

p +
  labs(title='Virus/Total reads ratio in infected cells', y="Virus/Total reads ratio", x = "Hours Post Infection") +
  theme_classic() +
  theme(legend.title = element_text(size = 20),
      legend.text = element_text(size = 15), 
      axis.text.x = element_text(face="bold", 
                                 size=14, angle = 45, hjust = 1),
      axis.text.y = element_text(face="bold", 
                                 size=14), 
      axis.title=element_text(size=14,face="bold") 
)

