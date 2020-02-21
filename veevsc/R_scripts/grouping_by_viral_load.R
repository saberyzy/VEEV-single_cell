# usage: mark cells based on viral load and generate the bar plot
# date updated: 02/19/2020
# author: Zhiyuan Yao
# other scripts needed : Preprocessing.R

library(Seurat)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)

# load data and mark cells based on viral load
VEEV_2 <- readRDS(file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")
metadata <- VEEV_2@meta.data
VEEV_2_nordata <- as.data.frame(x=as.matrix(x=VEEV_2@data))
metadata$viral_load_ridge <- 'uninfect'
high = row.names(subset(metadata, MOI != 0 & Virus_Total_ratio > 0.001))
low = row.names(subset(metadata, MOI != 0 & Virus_Total_ratio < 0.001))
metadata[high,'viral_load_ridge'] <- 'viral_load_high'
metadata[low,'viral_load_ridge'] <- 'viral_load_low'
VEEV_2@meta.data <- metadata

# calculate the percetage of different viral load groups in different time points
viruspercent <- data.frame('hour' = as.character(c(0, 0.5, 1.5, 4, 6, 12, 24)))
viruspercent$viral_load_low <- 0
viruspercent$viral_load_high <- 0
row.names(viruspercent) <- viruspercent$hour
VEEV_2 <- SetAllIdent(VEEV_2, id = "hours")
for (num in as.character(c(0.5, 1.5, 4, 6, 12, 24))){
  x = paste(num, "h", sep="")
  a <- length(row.names(subset(metadata, hours == x & viral_load_ridge == 'viral_load_low',)))
  b <- length(row.names(subset(metadata, hours == x & viral_load_ridge != 'uninfect',)))
  low <- (a/b) * 100
  high <- (1- a/b) * 100
  viruspercent[num, "viral_load_low"] <- low
  viruspercent[num, "viral_load_high"] <- high
}
viruspercent['0', "viral_load_low"] <- viruspercent['0.5', "viral_load_low"]
melted = melt(viruspercent, id.vars="hour")
colnames(melted) <- c('hour', 'viral_load', 'percentage')
melted$hour <- factor(melted$hour, levels = as.character(c(0, 0.5, 1.5, 4, 6, 12, 24)))

# plotting
p <- ggplot(melted) + geom_bar(aes(y = percentage, x = hour, fill = viral_load),                                stat="identity")
melted <- ddply(melted, .(hour), transform, pos = cumsum(percentage) - (0.5 * percentage))
p2 <- p + theme_classic() + scale_y_continuous(labels = dollar_format(suffix = "%", prefix = ""))
fill <- c("#5F9EA0", "#E1B378")
p3 <- p2 + scale_fill_manual(values=fill) + labs(x="Hours_post_infection", y="Percentage") +
  ggtitle("Viral load composition of infected cells (%)")
p4 <- p3 + 
  theme(legend.title = element_text(size = 18),
        legend.text = element_text(size = 12), 
        axis.text.x = element_text(face="bold", 
                                   size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(face="bold", 
                                   size=14), 
        axis.title=element_text(size=14,face="bold") 
  )
p4

saveRDS(VEEV_2, file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")


