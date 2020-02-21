# usage: calculate spearman rho and draw histogram plots
# date updated: 02/19/2020
# author: Zhiyuan Yao
# other scripts needed : Preprocessing.R

# load data
VEEV_2 <- readRDS(file = "~/workspace/VEEV/VEEV_Rproject/VEEV_2.rds")
metadata <- VEEV_2@meta.data
VEEV_2_nordata <- as.data.frame(x=as.matrix(x=VEEV_2@data))

# calculate spearman rho
infectcells <- row.names(subset(metadata, infect == 'Yes',))
VEEV_infect_ercc <-VEEV_2_nordata[, infectcells]
rownames_ercc <- row.names(VEEV_infect_ercc)
spearman_ercc <- data.frame(matrix(nrow=length(rownames_ercc), ncol=1), row.names = rownames_ercc)
colnames(spearman_ercc) <- 'spearman_cor'
VEEV_infect_ercc_T <- as.data.frame(t(VEEV_infect_ercc))

for (gene in rownames_ercc){
  targene <- VEEV_infect_ercc_T[,gene]
  testcor<- cor.test(targene, VEEV_infect_ercc_T$VEEVGFP, method="spearman", exact=FALSE)
  spearman_ercc[gene, 'spearman_cor'] <- testcor$estimate
}

spearman_ercc[is.na(spearman_ercc)] <- 0
spearman_ercc <- subset(spearman_ercc, spearman_cor != 0, )
spearman_ercc_var <- subset(spearman_ercc, row.names(spearman_ercc) != 'VEEVGFP',)

# plotting
p <- ggplot(spearman_ercc_var, aes(x=spearman_cor)) + 
  geom_histogram(aes(y=..density..),binwidth=.025, colour="black", fill="white") +  
  geom_density(alpha=.2, fill="#FF6666")
 
  p + labs(title="Histogram of VEEV-host corelation_ercc method", x = 'Correlation with VEEV',
           y = 'Density') +
    theme_classic() +
    theme(legend.title = element_text(size = 20),
          legend.text = element_text(size = 15), 
          axis.text.x = element_text(face="bold", 
                                     size=14),
          axis.text.y = element_text(face="bold", 
                                     size=14), 
          axis.title=element_text(size=14,face="bold"))
  
  
#save data
write.table(spearman_ercc, file = '/Users/yzhiyuan/workspace/VEEV/VEEV_Rproject/table/spearman_ercc.tsv', 
            quote = FALSE, sep = '\t')

