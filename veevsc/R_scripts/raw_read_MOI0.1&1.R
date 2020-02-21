# usage: plotting raw viral reads of cells at MOI = 1 or 0.1 at different timepoints
# date updated: 02/19/2020
# author: Zhiyuan Yao
# other scripts needed : Preprocessing.R

# load data
VEEV_2 <- readRDS(file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")
metadata <- VEEV_2@meta.data

# plotting
VEEV_2 <-SetAllIdent(VEEV_2, id = 'MOI')
moi <- c(0.1,1)
fig <-list()
for(x in moi){
print(x)
infectcells_moi <- row.names(subset(metadata, MOI == x, ))
print(length(infectcells_moi))
VEEV_moi <- SubsetData(object = VEEV_2, cells.use = infectcells_moi)
VEEV_moi <-SetAllIdent(VEEV_moi, id = 'hours')
subtitle <- paste('MOI = ', as.character(x), sep = '')

p <- VlnPlot(object = VEEV_moi, features.plot = c("VEEVGFP"), y.log = TRUE, 
             use.raw = TRUE, return.plotlist = TRUE)
p2 <- p + labs(x = 'Hours Post Infection (hpi)', y = 'Raw Viral Reads', 
         title = 'VEEV viral load at different time point', subtitle = subtitle) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 14, colour = 'black'), 
        axis.title.y = element_text(size = 14, colour = 'black'),
        plot.title = element_text(size = 14, colour = 'black')) 
fig[[as.character(x)]] <- p2
rm(VEEV_moi)
}

fig[['0.1']]
fig[['1']]


