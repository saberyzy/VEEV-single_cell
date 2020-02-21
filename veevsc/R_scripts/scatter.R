# usage: generate scatter plot for representative genes in spearman
# date updated: 02/19/2020
# author: Zhiyuan Yao
# other scripts needed: Preprocessing.R

library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)

# load the Rproject 
VEEV_2 <- readRDS(file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")
metadata <- VEEV_2@meta.data
VEEV_2_nordata <- as.data.frame(x=as.matrix(x=VEEV_2@data))
infectcells <- row.names(subset(metadata, infect == 'Yes',))

# normalize data with ERCC method with no log transfer in convenience of plotting 
VEEV_ercc <- VEEV_2@raw.data[, row.names(metadata)]
a  <- VEEV_ercc / (as.vector(metadata$ERCC_sum))
VEEV_2_ERCCnor <- a[row.names(VEEV_2_nordata),]
VEEV_2_nordata_infect <- VEEV_2_ERCCnor[,infectcells]

# plotting
genename <- c('TNFAIP3', 'TAF7', 'RND3')
spearman <- read.table('/Users/yzhiyuan/workspace/VEEV/VEEV_Rproject/table/spearman_ercc_1.tsv', 
                       sep='\t', row.names=1, header=T)
Scatter <- function(genename, spearman){
  plist <- list()
  s <- c(genename, 'VEEVGFP')
  df <- VEEV_2_nordata_infect[s,]
  df <- as.data.frame(t(df))
  df$hours <- metadata[infectcells, 'hours']
  print(colnames(df))
  for (x in genename){
    plist[[x]] <- ggplot(df, aes_string(x='VEEVGFP', y=x)) + 
      stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon') + 
      scale_fill_continuous(low="#d3d8de",high="#515457") +
      geom_point(aes(color = hours), alpha=0.3) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      labs(x = 'VEEV reads / ERCC reads', y = 'Gene reads / ERCC reads', title = x) +
      annotation_logticks() +
      theme_bw() +
      theme( 
        axis.text.x = element_text(face="bold", 
                                   size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(face="bold", 
                                   size=14),
        legend.position = 'none'
      )
    label <- paste('Rho = ',format(round(spearman[x,], 2), nsmall = 2), sep='')
    plist[[x]] <- plist[[x]] + annotate("text", x = 4, y = 25, label = label)
  }
  return(plist)
}

plist <- Scatter(genename, spearman)
quartz()
p <- grid.arrange(grobs = plist, ncol = length(plist))

ggsave('/Users/yzhiyuan/workspace/VEEV/VEEV_scripts/figures/figure_1E_grey.pdf', p, 
       width = 9, height = 3.22)

# try different KDE and contour setup
p <- ggplot(df, aes_string(x='VEEVGFP', y='TNFAIP3')) + 
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon') + 
  scale_fill_continuous(low="green",high="red") +
  geom_point(aes(color = hours), alpha=0.3) +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = 'VEEV reads / ERCC reads', y = 'TNFAIP3 reads / ERCC reads') +
  theme_bw() +
  theme( 
    axis.text.x = element_text(face="bold", 
                               size=14, angle = 45, hjust = 1),
    axis.text.y = element_text(face="bold", 
                               size=14),
    legend.position = 'none'
  )
p + annotate("text", x = 4, y = 25, label = "ρ = 0.25")








