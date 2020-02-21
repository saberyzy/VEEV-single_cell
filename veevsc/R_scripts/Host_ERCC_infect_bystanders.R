# usage: Host-ERCC ratio plots in infected and bystander cells
# date updated: 02/19/2020
# author: Zhiyuan Yao
# other scripts needed: Preprocessing.R


library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

# load data and Rprojects
VEEV_2 <- readRDS(file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")
metadata <- VEEV_2@meta.data
metadata$host_ERCC_ratio <- (metadata$Host_reads) / (metadata$ERCC_sum)

metadata_infect <- subset(metadata, infect == 'Yes', )
a <- metadata_infect[, c('hours', 'host_ERCC_ratio')]

# the plot for infected cells 
p <- ggplot(a, aes(hours,host_ERCC_ratio, fill = hours)) + 
  geom_violin(scale = 'width', trim = TRUE) +
  geom_boxplot(width=0.1, color="black", alpha=0.7, outlier.shape = NA, fill = "white")+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides = "l")

p +
  labs(title='Host/ERCC reads ratio in infected cells', y="Host/ERCC reads ratio", x = "Hours Post Infection") +
  theme_classic() +
  theme(legend.title = element_text(size = 20),
      legend.text = element_text(size = 15), 
      axis.text.x = element_text(face="bold", 
                                 size=14, angle = 45, hjust = 1),
      axis.text.y = element_text(face="bold", 
                                 size=14), 
      axis.title=element_text(size=14,face="bold") 
)


# the plot for bystanders
metadata_bystanders <- subset(metadata, infect == 'No' & MOI != 0, )
a <- metadata_bystanders[, c('hours', 'host_ERCC_ratio')]

p <- ggplot(a, aes(hours,host_ERCC_ratio, fill = hours)) + 
  geom_violin(scale = 'width', trim = TRUE) +
  geom_boxplot(width=0.1, color="black", alpha=0.7, outlier.shape = NA, fill = "white")+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides = "l")

p +
  labs(title='Host/ERCC reads ratio in bystander cells', y="Host/ERCC reads ratio", x = "Hours Post Infection") +
  theme_classic() +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 15), 
        axis.text.x = element_text(face="bold", 
                                   size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(face="bold", 
                                   size=14), 
        axis.title=element_text(size=14,face="bold") 
  )

#save
VEEV_2@meta.data <- metadata
saveRDS(VEEV_2, file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")