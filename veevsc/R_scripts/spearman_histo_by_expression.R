# usage: visualize spearman histogram with genes grouped by their expression
# date updated: 02/19/2020
# author: Zhiyuan Yao
# other scripts needed : spearman_histo.R

library(Seurat)
library(plyr)
library(dplyr)
library(ggplot2)

# load data
VEEV_2 <- readRDS(file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")
metadata <- VEEV_2@meta.data
VEEV_2_nordata <- as.data.frame(x=as.matrix(x=VEEV_2@data))

# select uninfect cells 
uninfect_cells <- row.names(subset(metadata, MOI == 0,))
VEEV_uninfect <- SubsetData(VEEV_2, cells = uninfect_cells)

# groups genes based on their expression
splitgroup <- function(object){
  namelist <- list()
  exp <- c(0, 2, 4, 6, 8)
  for (a in exp){
    if (a != 8){
    object <- FindVariableGenes(object = object, mean.function = ExpMean, dispersion.function = LogVMR, 
                                       x.low.cutoff = a, x.high.cutoff = a+2, y.cutoff = -Inf, y.high.cutoff = Inf, 
                                       do.plot = FALSE)
    }else{
      object <- FindVariableGenes(object = object, mean.function = ExpMean, dispersion.function = LogVMR, 
                                  x.low.cutoff = a, x.high.cutoff = Inf, y.cutoff = -Inf, y.high.cutoff = Inf, 
                                  do.plot = FALSE)
    }
    namelist[[as.character(a)]] <- object@var.genes
    
  }
  return(namelist)
}
name <- splitgroup(VEEV_uninfect)

# load spearman result and mark genes based on their expression
spearman_ercc <- read.table('/Users/yzhiyuan/workspace/VEEV/VEEV_Rproject/table/spearman_ercc_1.tsv', sep='\t', 
                             row.names=1, header=T)
spearman_ercc_exp <- subset(spearman_ercc, spearman_cor != 0, )

spearman_ercc[name[['0']],'expression'] <- 'very low'
spearman_ercc[name[['2']],'expression'] <- 'low'
spearman_ercc[name[['4']],'expression'] <- 'median'
spearman_ercc[name[['6']],'expression'] <- 'high'
spearman_ercc[name[['8']],'expression'] <- 'very high'

spearman_ercc_exp <- spearman_ercc[unlist(name, use.names=FALSE),]
spearman_ercc_exp <- spearman_ercc_exp[complete.cases(spearman_ercc_exp), ]
spearman_ercc_exp <- subset(spearman_ercc_exp, spearman_cor != 0, )
spearman_ercc_exp$expression <- factor(spearman_ercc_exp$expression, 
                                   levels = c("very low", "low", "median", "high", "very high"))
aver <- ddply(spearman_ercc_exp, "expression", summarise, median=median(spearman_cor))

# plotting
p <- ggplot(spearman_ercc_exp, aes(x=spearman_cor, fill=expression)) + 
  geom_density(alpha=.3) +
  geom_vline(data=aver, aes(xintercept=median,  colour=expression),
             linetype="dashed", size=1) 

p + labs(title="Histogram of VEEV-host corelation", x = 'Correlation with VEEV',
     y = 'Density') +
  theme_classic() +
  theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 15), 
        axis.text.x = element_text(face="bold", 
                                   size=14),
        axis.text.y = element_text(face="bold", 
                                   size=14), 
        axis.title=element_text(size=14,face="bold"))



