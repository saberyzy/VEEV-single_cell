# usage: calculate and plot infection rates
# date updated: 02/19/2020
# author: Zhiyuan Yao
# other scripts needed : Preprocessing.R

library(Seurat)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)

# Load Rprojects
VEEV_2 <- readRDS(file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")
metadata <- VEEV_2@meta.data
VEEV_2_nordata <- as.data.frame(x=as.matrix(x=VEEV_2@data))

# calculate the infection rate
MOIfraction <- list()
for (moi in c(0, 0.1, 1)){
  MOIfraction[[as.character(moi)]] <- numeric()
  for(timepoint in c('0.5h', '1.5h', '4h', '6h', '12h', '24h')){
    infect <- length(rownames(subset(metadata, infect == 'Yes' & MOI == moi & hours == timepoint, )))
    uninfect <- length(rownames(subset(metadata, infect == 'No' & MOI == moi & hours == timepoint, )))
    Frac <- infect/(infect + uninfect)
    MOIfraction[[as.character(moi)]] <- c(MOIfraction[[as.character(moi)]], Frac)
    print(MOIfraction[[as.character(moi)]])
  }
}
MOIfractable <- as.data.frame(MOIfraction)
row.names(MOIfractable) <- c('0.5h', '1.5h', '4h', '6h', '12h', '24h')
MOIfractable$timepoint <- row.names(MOIfractable)
colnames(MOIfractable) <-c('0', '0.1', '1', 'timepoint')
MOIfractable$timepoint <- ordered(MOIfractable$timepoint, levels = MOIfractable$timepoint[c(1, 2, 3, 4, 5, 6)])
levels(MOIfractable$timepoint)

#plotting
melted = melt(MOIfractable, id.vars="timepoint")
colnames(melted) <- c('Infect_hours', 'MOI', 'Infect_fraction')
p1 <- ggplot(melted, aes(x = Infect_hours,  y=Infect_fraction, group=MOI, shape=MOI, color=MOI)) + 
  geom_point() + 
  geom_line() +
  scale_colour_manual(values=c('#999999','#E69F00', '#56B4E9'))
p1 + 
  labs(title="VEEV Infection", y="Infection rate", x = "Hours post infection") +
  theme_classic() +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 15), 
        axis.text.x = element_text(face="bold", 
                                   size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(face="bold", 
                                   size=14), 
        axis.title=element_text(size=14,face="bold") 
  )


  




