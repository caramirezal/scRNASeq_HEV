## Mapping gene clusters signatures obtain from a infected human cell line

## Dependencies
library(Seurat)
library(dplyr)
library(ggplot2)
library(AUCell)
library(readr)
library(readxl)
library(viridis)
library(reshape2)
library(ggside)
library(tidyquant)
library(ggExtra)
library(ggpubr)
library(cowplot)
library(enrichR)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)


##--------------------------------------------------------------------------
## Initial settings

set.seed(333)

rerun <- FALSE

## Working directory
path2project <- '/media/ag-cherrmann/cramirez/scRNASeq_HEV'
setwd(path2project)


## Loading general settings
source(paste0(path2project, '/scripts/R/settings.R'))


path2figures <- paste0(path2project, '/figures/04_gene_clusters')
if ( ! dir.exists(path2figures)) {
    dir.create(path2figures)
}

path2analysis <- paste0(path2project, '/analysis/04_gene_clusters')
if ( ! dir.exists(path2analysis)) {
    dir.create(path2analysis)
}


## Reading seurat objects
seurat_2 <- readRDS('/media/ag-cherrmann/projects/24_scRNA_Loan/analysis/seuratIntegrated_runsJan2024_filtprocd.rds')
seurat_2



## Loading clusterWT genes analysis
path2clusterGenes <- '/media/ag-cherrmann/lcosta/Interferon/data/clusterWTgenes_CmeansMicroarray2023-10-19.rds'
clusterGenes <- readRDS(path2clusterGenes)
names(clusterGenes)
head(clusterGenes$C1)


## Saving genes into xlsx
writexl::write_xlsx(clusterGenes, paste0(path2analysis, "/gene_clusters_microarray_Anna.xlsx"))




if ( rerun == TRUE ) {
    counts_matrix <- seurat_2@assays$RNA@layers$counts
    rownames(counts_matrix) <- rownames(seurat_2)
    colnames(counts_matrix) <- colnames(seurat_2)
    dim(counts_matrix)
    counts_matrix[1:5, 1:5]

    ## Calculating the scores
    cells_rankings <- AUCell_buildRankings(counts_matrix, plotStats = FALSE)
    cells_AUC <- AUCell_calcAUC(signatures_list, cells_rankings, 
                                aucMaxRank=(nrow(cells_rankings)*0.05))
    auc.mtx <- getAUC(cells_AUC)
    auc.df <- as.data.frame(t(auc.mtx))
    auc.df$"barcode" <- rownames(auc.df)
    write_tsv(auc.df, file = paste0(path2analysis, "/gene_clusters_aucell.tsv.gz"))
}


## Check barcodes order
all(colnames(seurat_2) == rownames(auc.df))
seurat_2@meta.data <- cbind(seurat_2@meta.data, auc.df)




## function to plot features
plot_umap <- function(
        signature=signature_name,
        point_size=0.5,
        quant = 0.95
        ){ 
        barcodes <- seurat_2@meta.data[signature_name] %>%
            filter(!!as.symbol(signature_name) < 
                        quantile(!!as.symbol(signature_name), quant) ) %>%
                        rownames()
        umap <- FeaturePlot(
                    seurat_2, 
                    cells = barcodes,
                    reduction = "umap.integrated.cca", 
                    features = signature,
                    pt.size = point_size,
                    order=TRUE) +
                theme_bw() +
                theme(panel.border = element_blank()) +
                scale_colour_viridis() +
                labs(title=signature)

        return(umap)
}


signature_names <- paste0("C", 1:4)
umap_list <- lapply(
    signature_names, 
    function(signature_name){
        plot_umap(signature = signature_name, point_size = 0.1)
})



## @plot umap
pdf(paste0(path2figures, "/umaps_gene_clusters.pdf"),
    height = 10, width = 10)
gridExtra::grid.arrange(grobs=umap_list, ncol=2)
dev.off()



##-------------------------------------------------------------
## boxplots signatures by condition

signatures_df <- FetchData(
    seurat_2, 
    vars = c(paste0("C", 1:4), 
            "umapintegratedcca_1",
            "umapintegratedcca_2",
            "orig.ident", 
            "infected"),
    slot="data",
)
head(signatures_df)


pdf(paste0(path2figures, "/boxplots_gene_clusters_signatures.pdf"),
    height = 8, width = 8)
boxplots <- signatures_df  %>%
             select(-umapintegratedcca_1, -umapintegratedcca_2) %>%
             melt() %>%
                mutate(orig.ident=factor(orig.ident, 
                            level=sort(unique(signatures_df$"orig.ident")))) %>%
                ggboxplot(x="orig.ident",
                          y="value",
                          fill="infected") +
                        facet_wrap(~variable, scale="free") +
                        scale_fill_manual(values=colors_infection) +
                        theme_bw() +
                        stat_compare_means(aes(group = infected),
                                           label = "p.signif") +
                        labs(x="", y="AUC Score") 
boxplots
dev.off()