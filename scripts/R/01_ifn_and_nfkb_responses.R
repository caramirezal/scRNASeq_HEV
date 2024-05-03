## Analysis of the transcriptomic response to HEV infection in human liver cell lines
## for early and late time points and virus carrying mutations for ORF genes



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


path2figures <- paste0(path2project, '/figures/01_ifn_and_nfkb_responses')
if ( ! dir.exists(path2figures)) {
    dir.create(path2figures)
}

path2analysis <- paste0(path2project, '/analysis/01_ifn_and_nfkb_responses')
if ( ! dir.exists(path2analysis)) {
    dir.create(path2analysis)
}


## Reading seurat objects
#seurat_1 <- readRDS('/media/ag-cherrmann/projects/24_scRNA_Loan/data/run_Feb2023/HEVmergedConds_normalised_all_annot.rds')
seurat_2 <- readRDS('/media/ag-cherrmann/projects/24_scRNA_Loan/analysis/seuratIntegrated_runsJan2024_filtprocd.rds')





## Inspecting metadata
#head(seurat_1)
head(seurat_2)


##--------------------------------------------------------
## Evaluation of IFN and NFkB signatures

## Extracting counts matrix

if ( rerun == TRUE ) {
    counts_matrix <- seurat_2@assays$RNA@layers$counts
    rownames(counts_matrix) <- rownames(seurat_2)
    colnames(counts_matrix) <- colnames(seurat_2)
    dim(counts_matrix)
    counts_matrix[1:5, 1:5]
}



## Reading the signatures
ifn_signature_df <- readxl::read_xlsx(
    paste0(path2project, 
    '/data/signatures/240412_ISG-List_Schoggins_FINAL.xlsx'),
    col_names = FALSE
)
nfkb_signature_df <- readxl::read_xlsx(
    paste0(path2project, 
    '/data/signatures/240412_NF-kB-Gene-List_FINAL.xlsx'),
    col_names = FALSE
)
signatures_list <- list(ifn_signature=ifn_signature_df$'...1',
                        nfkb_signature=nfkb_signature_df$'...1')


if ( rerun == TRUE ) {
    ## Calculating the scores
    cells_rankings <- AUCell_buildRankings(counts_matrix, plotStats = FALSE)
    cells_AUC <- AUCell_calcAUC(signatures_list, cells_rankings, 
                                aucMaxRank=(nrow(cells_rankings)*0.05))
    auc.mtx <- getAUC(cells_AUC)
    auc.df <- as.data.frame(t(auc.mtx))
    auc.df$"barcode" <- rownames(auc.df)
    write_tsv(auc.df, file = paste0(path2analysis, "/ifn_nfkb_aucell.tsv.gz"))
}


auc.df <- read_tsv(paste0(path2analysis, "/ifn_nfkb_aucell.tsv.gz"))
auc.df <- as.data.frame(auc.df)
rownames(auc.df) <- auc.df$"barcode" 
head(auc.df)
seurat_2 <- AddMetaData(seurat_2, auc.df)



## extracting embeddings and metadata
gex_df <- FetchData(
    seurat_2, 
    vars = c(colnames(head(seurat_2)))
)
umap_df <- Embeddings(seurat_2, 
reduction='umap.integrated.cca')
gex_df <- cbind(gex_df, umap_df)

## Inspecting samples
seurat_2$orig.ident %>% table() 




##--------------------------------------------------------------------------
## Plotting UMAP of the integrated/unintegrated samples
pdf(paste0(path2figures, '/umap_integrated_samples_seq2.pdf'),
    width=10, height=4.5)
umap_1 <- DimPlot(seurat_2, 
        group.by = 'orig.ident',
        cols=colors_seq2_samples) +
        theme_bw() +
        theme(panel.border = element_blank()) +
        ggtitle('')
umap_2 <- DimPlot(seurat_2, 
        reduction = "umap.integrated.cca", 
        group.by = 'orig.ident',
        cols=colors_seq2_samples) +
        theme_bw() +
        theme(panel.border = element_blank()) +
        ggtitle('')
umap_1 + umap_2
dev.off()





##--------------------------------------------------------------
## QC metrics
umap_nCount <- FeaturePlot(seurat_2, 
                           features='nCount_RNA',
                           order=TRUE) +
        theme_bw() +
        theme(panel.border = element_blank()) +
        ggtitle('') +
        scale_color_viridis()

umap_nFeature <- FeaturePlot(seurat_2, 
                           features='nFeature_RNA',
                           order=TRUE) +
        theme_bw() +
        theme(panel.border = element_blank()) +
        ggtitle('') +
        scale_color_viridis()



pdf(paste0(path2figures, '/umap_qc_metrics.pdf'),
     width = 4.5, height=4.5)
umap_nCount + umap_nFeature
dev.off()


jpeg(paste0(path2figures, '/umap_qc_metrics.jpeg'))
umap_nCount + umap_nFeature
dev.off()







## Plotting UMAP of the integrated samples
pdf(paste0(path2figures, '/umap_ifn_nfkb_signatures.pdf'),
    width=10, height=4.5)
umap_ifn <- FeaturePlot(subset(seurat_2,
                               ifn_signature < quantile(seurat_2$"ifn_signature", 0.95)), 
                      features = 'ifn_signature',
                      order=TRUE) +
        theme_bw() +
        scale_colour_viridis() +
        theme(panel.border = element_blank()) +
        ggtitle('IFN')
umap_nfkb <- FeaturePlot(subset(seurat_2,
                               ifn_signature <  quantile(seurat_2$"ifn_signature", 0.95)), 
                         features = 'nfkb_signature',
                         order=TRUE) +
        theme_bw() +
        scale_colour_viridis() +
        theme(panel.border = element_blank()) +
        ggtitle('NFkB')
umap_ifn | umap_nfkb
dev.off()



## Plotting UMAP of the integrated samples
pdf(paste0(path2figures, '/umap_ifn_nfkb_signatures_integrated.pdf'),
    width=10, height=4.5)
umap_ifn <- FeaturePlot(subset(seurat_2,
                               ifn_signature < quantile(seurat_2$"ifn_signature", 0.95)), 
                      features = 'ifn_signature',
                      reduction = "umap.integrated.cca", 
                      order=TRUE) +
        theme_bw() +
        scale_colour_viridis() +
        theme(panel.border = element_blank()) +
        ggtitle('IFN')
umap_nfkb <- FeaturePlot(subset(seurat_2,
                               ifn_signature <  quantile(seurat_2$"ifn_signature", 0.95)), 
                         features = 'nfkb_signature',
                         reduction = "umap.integrated.cca", 
                         order=TRUE) +
        theme_bw() +
        scale_colour_viridis() +
        theme(panel.border = element_blank()) +
        ggtitle('NFkB')
umap_ifn | umap_nfkb
dev.off()



jpeg(paste0(path2figures, '/umap_ifn_nfkb_signatures.jpg'))
umap_ifn | umap_nfkb
dev.off()




##--------------------------------------------------------------------------
## Viral gene expression orfs

jpeg(paste0(path2figures, '/boxplot_orfs.jpg'))
boxplot_orfs <- gex_df %>%
    select(ORF1, ORF2, ORF3, orig.ident) %>%
    melt() %>%
    ggplot(aes(x=orig.ident,
               y=log(value+1),
               fill=orig.ident)) +
        geom_boxplot() +
        theme_bw() +
        facet_wrap(~variable, 
                   scale='free') +
        scale_fill_manual(values=colors_seq2_samples) +
        theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank())
boxplot_orfs
dev.off()




pdf(paste0(path2figures, '/boxplot_orfs.pdf'),
    height = 3.5, width = 7)
boxplot_orfs <- gex_df %>%
    select(ORF1, ORF2, ORF3, orig.ident) %>%
    melt() %>%
    ggplot(aes(x=orig.ident,
               y=log(value+1),
               fill=orig.ident)) +
        geom_boxplot() +
        theme_bw() +
        facet_wrap(~variable, 
                   scale='free') +
        scale_fill_manual(values=colors_seq2_samples) +
        theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank()) +
        labs(fill='')
boxplot_orfs
dev.off()




##------------------------------------------------------------------
## Distribution of ORF1/2
jpeg(paste0(path2figures, '/scatterplot_orf_distribution.jpeg'))
scatterplot_orf_distribution <- gex_df %>%
    ggplot(aes(x=ORF1+1, 
               y=ORF2+1, 
               colour=infected)) +
            geom_point(alpha=0.5) +
            theme_classic() + 
            scale_x_continuous(trans='log10') +
            scale_y_continuous(trans='log10') +
            facet_wrap(~orig.ident)
scatterplot_orf_distribution
dev.off()



pdf(paste0(path2figures, 
             '/scatterplot_orf_distribution.pdf'),
    height=6, width = 6.5)
scatterplot_orf_distribution 
dev.off()






##--------------------------------------------------------------------------
## Bystander cell definition
## Distribution of INF and NFkB scores in cells

jpeg(paste0(path2figures, '/scatterplot_orf2_vs_ifn.jpg'))
scatterplot_orf2_ifn <- gex_df %>% 
    ggplot(aes(x=ORF2+1,
               y=ifn_signature,
               colour=orig.ident)) +
               geom_point(alpha=0.75) +
               theme_bw() +
            scale_x_continuous(trans='log10') +
            scale_colour_manual(values=colors_seq2_samples) +
            labs(colour='',
                 y='IFN Signature')
scatterplot_orf2_ifn     
dev.off()


pdf(paste0(path2figures, '/scatterplot_orf2_vs_ifn.pdf'),
    width = 8, height = 5.5)
scatterplot_orf2_ifn
dev.off()



## Exploration of the definition of IFN active cells

## Histogram of IFN expression scores
head(gex_df) 
thresholds <- c(0.04, 0.06, 0.08, 0.1, 0.12, 0.14)
jpeg(paste0(path2figures, "/histogram_ifn_scores.jpg"))
gex_df %>%
    ggplot(aes(x=ifn_signature,
               fill=orig.ident)) +
        geom_density(alpha=0.3) +
        geom_vline(xintercept = thresholds,
                   linetype="dashed",
                   colour="red") +
        facet_wrap(~orig.ident) +
        theme_bw() +
        scale_fill_manual(values=colors_seq2_samples)
dev.off()



th <- thresholds[1]
seurat_2$"ifn_level" <- ifelse(seurat_2$"ifn_signature">th, 
                                "high", "low")


## Adding unintegrated umap info
umap_unintegrated <- seurat_2@reductions$umap.unintegrat@cell.embeddings
all(rownames(umap_unintegrated) == rownames(gex_df))
gex_df <- cbind(gex_df, umap_unintegrated)
gex_df$"ifn_level" <- seurat_2$"ifn_level"


plot_umap_threshold <- function(th){
    gex_df$"ifn_level" <- ifelse(gex_df$"ifn_signature">th, 
                                "high", "low")
    gex_df %>%
        filter(ifn_level=="low") %>%
        ggplot(aes(x=umapintegratedcca_1, 
                y=umapintegratedcca_2)) +
                geom_point(size=0.1,
                            colour="gray90") +
                    geom_point(data=filter(gex_df, 
                                           ifn_level=="high"),
                                mapping=aes(x=umapintegratedcca_1, 
                                            y=umapintegratedcca_2),
                                colour="red",
                                size=0.1) +
                theme_bw() +
                theme(panel.border = element_blank()) +
                labs(title=paste0("Threshold = ", th)) 
}
names(thresholds) <- paste0("threshold=", thresholds)
umap_threshold_list <- lapply(thresholds, plot_umap_threshold)



jpeg(paste0(path2figures, "/umap_threshold.jpg"))
gridExtra::grid.arrange(grobs=umap_threshold_list, ncol=3)
dev.off()


isgs <- c("LGALS3", "ISG15", "B2M", "RPL22", "STAT1",  "IRF1")
isgs_gex <- FetchData(seurat_2, vars=isgs)
#isgs_gex <- FetchData(seurat_2, vars="RPL22")
gex_df <- cbind(gex_df, isgs_gex)


names(isgs) <- isgs
plot_umap_markers <- function(marker){
    gex_df <- arrange(gex_df, !!as.name(marker))
    gex_df %>%
        ggplot(aes(x=umapintegratedcca_1, 
                   y=umapintegratedcca_2,
                   colour=!!as.name(marker))) +
                geom_point(size=0.01) +
                theme_bw() +
                theme(panel.border = element_blank()) +
                labs(title=marker) +
                scale_colour_viridis(option="magma")
}
umap_markers_list <- lapply(isgs, plot_umap_markers)


jpeg(paste0(path2figures, "/umap_isgs.jpg"))
gridExtra::grid.arrange(grobs=umap_markers_list, ncol=3)
dev.off()


pdf(paste0(path2figures, "/umap_isgs.pdf"),
    height = 8, width = 10)
gridExtra::grid.arrange(grobs=umap_markers_list, ncol=3)
dev.off()



jpeg(paste0(path2figures, "/umap_negative_cells.jpg"))
gex_df <- arrange(gex_df, desc(ifn_signature))
umap_uninfected <- gex_df %>%
        ggplot(aes(x=umapunintegrated_1, 
                y=umapunintegrated_2)) +
                geom_point(size=0.1,
                            colour="gray90") +
                    geom_point(data=filter(gex_df, 
                                           ORF2<1 &
                                            ! grepl("Mock", orig.ident) &
                                           ifn_signature < 0.0671),
                                mapping=aes(x=umapunintegrated_1, 
                                            y=umapunintegrated_2,
                                            colour=ifn_signature),
                                size=0.1) +
                    scale_color_viridis(option="magma") +
                theme_bw() +
                theme(panel.border = element_blank()) 
umap_uninfected
dev.off()



pdf(paste0(path2figures, "/umap_negative_cells.pdf"))
umap_uninfected
dev.off()


gex_df$"ifn_level" <- ifelse(gex_df$"ifn_signature">0.08, 
                                "high", "low")



plot_umap_heterogeneity <- function(
    seurat_subset,
    sample,
    point_size=1,
    min.dist = 0.0001,
    n.neighbors = 600,
    signature = "ifn_signature",
    resolution = 0.1,
    infection_status = "all" 
){
    seurat_subset <- subset(seurat_subset[signatures_list[signature][[1]], ], 
                                orig.ident == sample )
    cts <- seurat_subset@assays$RNA@layers$counts
    colnames(cts) <- colnames(seurat_subset)

    seurat <- CreateSeuratObject(
                counts = cts,
                assay = "RNA", 
                min.cells=0, 
                min.features=0, project="HEV", 
                meta.data = seurat_subset@meta.data
    )
    seurat <- NormalizeData(seurat)
    seurat <- FindVariableFeatures(seurat)
    seurat <- ScaleData(seurat)
    seurat <- RunPCA(seurat)
    seurat <- RunUMAP(seurat, 
                    dims = 1:30, 
                    reduction = "pca", 
                    reduction.name = "umap.unintegrated", 
                    min.dist = min.dist,
                    n.neighbors = n.neighbors)
    seurat <- FindNeighbors(seurat)
    seurat <- FindClusters(seurat, resolution=resolution)

    umap_coord_df <- FetchData(seurat, 
                        vars = c("umapunintegrated_1",
                                    "umapunintegrated_2",
                                    signature,
                                    "infected",
                                    "ORF2",
                                    "seurat_clusters"))
    umap_coord_df$"barcode" <- colnames(seurat)

    ## ifn signature
    if ( grepl("Mock", sample)) { quant <- 1 } else {  quant <-  quantile(as.vector(umap_coord_df[,signature]), 0.99) }
    umap_ifn <- umap_coord_df %>%
        filter(!!as.name(signature) < quant )  %>%
        arrange(!!as.name(signature)) %>%
        ggplot(aes(x=umapunintegrated_1, 
                y=umapunintegrated_2,
                colour=!!as.name(signature))) +
                geom_point(size=point_size) +
                theme_classic() +
                scale_colour_viridis() +
                labs(title=paste0(sample, " - ", signature))

    ## orf2
    if (infection_status == "all" ) {
        if ( grepl("Mock", sample)) { quant <- 1 } else {  quant <-  quantile(umap_coord_df$ORF2, 0.95) }
        umap_orf <- umap_coord_df %>%
            filter(ORF2 < quant ) %>%
            arrange(ORF2) %>%
            ggplot(aes(x=umapunintegrated_1, 
                    y=umapunintegrated_2,
                    colour=ORF2)) +
                    geom_point(size=point_size) +
                    theme_classic() +
                    scale_colour_viridis(option="magma") +
                    labs(title=paste0(sample, " - ", "ORF2"))
    }
    ## clusters
    colors_clusters <- RColorBrewer::brewer.pal(name="Set2", 
                                                n=8)[1:length(unique(umap_coord_df$seurat_clusters))]
    names(colors_clusters) <- as.character( 0:(length(colors_clusters)-1) )

    if ( sample == 'WTEl' ) { colors_clusters[c(1,2)] = colors_clusters[c(2,1)] }
    if ( grepl('Mock', sample) ) { colors_clusters[1] = "gray40" }

    umap_clusters <- umap_coord_df %>%
        ggplot(aes(x=umapunintegrated_1, 
                   y=umapunintegrated_2,
                   colour=seurat_clusters)) +
                geom_point(size=point_size) +
                scale_colour_manual(values=colors_clusters) +
                theme_classic() +
                labs(title=paste0(sample, " - ", "Clusters"))
    ## Saving results
    write_tsv(umap_coord_df, paste0(path2analysis, 
                              '/umap_coord_', sample,
                              "_", signature, "_", infection_status,
                              '.tsv.gz'))

    if ( infection_status == "all") {
        res = list(umap_ifn=umap_ifn, umap_orf=umap_orf,
                umap_clusters=umap_clusters) }
    if ( infection_status == "negative") {
        res = list(umap_ifn=umap_ifn, 
                umap_clusters=umap_clusters) }
    res
}





##-------------------------------------------------------------------------------------------
## All cells 

## IFN signature

## Clustering
samples <- unique(seurat_2$orig.ident)
#samples <- "WTEl"
umaps_ifn_heterogeneity_list <- lapply(
    sort(samples), function(samp){
            plot_umap_heterogeneity(seurat_subset = seurat_2,
                                    sample = samp,
                                    signature = "ifn_signature",
                                    infection_status = "all")
    }
)
umaps_ifn_heterogeneity_unlisted <- unlist(umaps_ifn_heterogeneity_list, recursive = FALSE)
names(umaps_ifn_heterogeneity_list) <- samples
jpeg(paste0(path2figures, "/umap_heterogeneity_only_ifn_genes_all.jpg"))
gridExtra::grid.arrange(
    grobs=umaps_ifn_heterogeneity_unlisted, 
    ncol=6
)
dev.off()


pdf(paste0(path2figures, "/umap_heterogeneity_only_ifn_genes_all.pdf"),
    width=18, height = 10)
gridExtra::grid.arrange(
    grobs=umaps_ifn_heterogeneity_unlisted, 
    ncol=6
)
dev.off()



##--------------------------------------------------------------
## Negative cells

## IFN signature

## Clustering
seurat_negative <- subset(seurat_2, infected == FALSE)
samples <- unique(seurat_2$orig.ident)
#samples <- "WTEl"
umaps_ifn_heterogeneity_list <- lapply(
    sort(samples), function(samp){
            plot_umap_heterogeneity(seurat_subset = seurat_negative,
                                    sample = samp,
                                    signature = "ifn_signature",
                                    infection_status = "negative")
    }
)
umaps_ifn_heterogeneity_unlisted <- unlist(umaps_ifn_heterogeneity_list, recursive = FALSE)
names(umaps_ifn_heterogeneity_list) <- samples
jpeg(paste0(path2figures, "/umap_heterogeneity_only_ifn_genes_negative.jpg"))
gridExtra::grid.arrange(
    grobs=umaps_ifn_heterogeneity_unlisted, 
    ncol=6
)
dev.off()


pdf(paste0(path2figures, "/umap_heterogeneity_only_ifn_genes_negative.pdf"),
    width=14, height = 8.5)
gridExtra::grid.arrange(
    grobs=umaps_ifn_heterogeneity_unlisted, 
    ncol=4
)
dev.off()




##--------------------------------------------------------------



## NFkB signature

## Clustering
samples <- unique(seurat_2$orig.ident)
umaps_nfkb_heterogeneity_list <- lapply(
    sort(samples), function(samp){
            plot_umap_heterogeneity(seurat_subset = seurat_2,
                                    sample = samp,
                                    signature = "nfkb_signature",
                                    infection_status = "all")
    }
)
umaps_nfkb_heterogeneity_unlisted <- unlist(umaps_nfkb_heterogeneity_list, recursive = FALSE)
names(umaps_nfkb_heterogeneity_list) <- samples
jpeg(paste0(path2figures, "/umap_heterogeneity_only_nfkb_genes_all.jpg"))
gridExtra::grid.arrange(
    grobs=umaps_nfkb_heterogeneity_unlisted, 
    ncol=6
)
dev.off()


pdf(paste0(path2figures, "/umap_heterogeneity_only_nfkb_genes_all.pdf"),
    width=18, height = 10)
gridExtra::grid.arrange(
    grobs=umaps_nfkb_heterogeneity_unlisted, 
    ncol=6
)
dev.off()



##-------------------------------------------------------------------


## Clustering
seurat_negative <- subset(seurat_2, infected == FALSE)
samples <- unique(seurat_2$orig.ident)
umaps_ifn_heterogeneity_list <- lapply(
    sort(samples), function(samp){
            plot_umap_heterogeneity(seurat_subset = seurat_negative,
                                    sample = samp,
                                    signature = "nfkb_signature",
                                    infection_status = "negative")
    }
)
umaps_ifn_heterogeneity_unlisted <- unlist(umaps_ifn_heterogeneity_list, recursive = FALSE)
names(umaps_ifn_heterogeneity_list) <- samples
jpeg(paste0(path2figures, "/umap_heterogeneity_only_nfkb_genes_negative.jpg"))
gridExtra::grid.arrange(
    grobs=umaps_ifn_heterogeneity_unlisted, 
    ncol=6
)
dev.off()


pdf(paste0(path2figures, "/umap_heterogeneity_only_nfkb_genes_negative.pdf"),
    width=14, height = 8.5)
gridExtra::grid.arrange(
    grobs=umaps_ifn_heterogeneity_unlisted, 
    ncol=4
)
dev.off()






##-------------------------------------------------------------------------------------------



signature <- 'nfkb_signature'
sample = "WTEl"
seurat_negative = subset(seurat_2, infected == FALSE)
plot_heatmap = function(sample = "WTEl",
                         signature = "ifn_signature",
                         seurat_subset,
                         infection_status = 'negative' ){
    signature_df = FetchData(subset(seurat_subset, orig.ident==sample), 
                                    vars=c(signatures_list[signature][[1]], 
                                            "ORF2", "ifn_signature",
                                            "nfkb_signature"), 
                                    layer="data")
    signature_df = mutate(signature_df, ORF2_log=log10(ORF2+1))
    cluster_anns = read_tsv(paste0(path2analysis, 
                                '/umap_coord_', sample,
                                "_", signature, 
                                '_', infection_status,
                                '.tsv.gz'))
    signature_df$"seurat_clusters" <- as.character(cluster_anns$"seurat_clusters")
    mtx = select(signature_df, 
                -ORF2, 
                -ORF2_log, 
                -ifn_signature, 
                -seurat_clusters, 
                -nfkb_signature) %>%
                    as.matrix() %>%
                    t()
    hvgs = apply(mtx, 1, var) %>%
                sort(decreasing = TRUE) %>%
                head(50) %>%
                names()
    mtx = mtx[hvgs, ]
    mdata = select(signature_df, 
                    ORF2, 
                    ORF2_log,
                    seurat_clusters, (!!as.name(signature)))
    min <- min(pull(mdata, !!as.name(signature)))
        max = max(pull(mdata, !!as.name(signature)))
        col_fun_ifn = circlize::colorRamp2(breaks=c(min,
                                                    ( max - min ) / 2,
                                                    max), 
                                            c("blue",
                                            "yellow",
                                            "red"))
    
    colors_clusters <- RColorBrewer::brewer.pal(name="Set2", n=8)[1:length(unique(mdata$"seurat_clusters"))]
    names(colors_clusters) <- unique(mdata$seurat_clusters) 
    if ( sample == 'WTEl') { 
        colors_clusters[c(1,2)] = colors_clusters[c(2,1)] 
    }
    if ( grepl('Mock', sample) ) { 
        colors_clusters[1] = "gray40" 
    }
    if ( infection_status == 'negative' ) {
        col_ann <- HeatmapAnnotation(Signature=mdata[signature][[1]],
                                seurat_clusters=mdata$seurat_clusters,
                                col = list(Signature=col_fun_ifn,
                                            seurat_clusters=colors_clusters
                                            ))
    }
    if ( infection_status == 'all' ) {
        min <- min(pull(mdata, !!as.name(signature)))
        max <- max(pull(mdata, !!as.name(signature)))
        col_fun_ifn = circlize::colorRamp2(breaks=c(min,
                                                    ( max - min ) / 2,
                                                    max), 
                                            c("blue",
                                            "yellow",
                                            "red"))
        min <- min(pull(mdata, ORF2_log))
        max <- max(pull(mdata, ORF2_log))
        col_fun_orf2 = circlize::colorRamp2(breaks=c(min,
                                                    ( max - min ) / 2,
                                                    max), 
                                            c("blue",
                                            "yellow",
                                            "red"))
        col_ann = HeatmapAnnotation(Signature=mdata[signature][[1]], 
                        ORF2_log=mdata$ORF2_log,
                        seurat_clusters=mdata$seurat_clusters,
                        col = list(Signature=col_fun_ifn,
                                    ORF2_log=col_fun_orf2,
                                    seurat_clusters=colors_clusters))
    } 
    Heatmap(mtx, top_annotation=col_ann, show_column_names=FALSE)
}

## ORF2 negative cells

samples <- unique(seurat_2$orig.ident) %>% sort()
signatures <- names(signatures_list)
for ( sign in  signatures ) {
    for ( samp in samples) {
        pdf(paste0(path2figures, "/heatmap_", 
                   samp, '_', sign, "_negative.pdf"),
             width=10, height = 16)
        hm <- plot_heatmap(sample = samp, 
                           signature = sign, 
                           seurat_subset = seurat_negative)
        draw(hm)
        dev.off()
    }
}


samples <- unique(seurat_2$orig.ident) %>% sort()
samples
signatures <- names(signatures_list)
for ( sign in  signatures ) {
    for ( samp in samples) {
        pdf(paste0(path2figures, "/heatmap_", 
                   samp, '_', sign, "_all.pdf"),
             width=10, height = 16)
        hm <- plot_heatmap(sample = samp, 
                           signature = sign, 
                           seurat_subset = seurat_2,
                           infection_status = 'all')
        draw(hm)
        dev.off()
    }
}




##------------------------------------------------------------------------
## Only using negative cells


## IFN signature

## Clustering
seurat_infected <- subset(seurat_2, infected == FALSE)
samples <- unique(seurat_2$orig.ident)
#samples <- "WTEl"
umaps_ifn_heterogeneity_list <- lapply(
    sort(samples), function(samp){
            plot_umap_heterogeneity(seurat_subset = seurat_infected,
                                    sample = samp,
                                    signature = "ifn_signature",
                                    infection_status = "negative")
    }
)
umaps_ifn_heterogeneity_unlisted <- unlist(umaps_ifn_heterogeneity_list, recursive = FALSE)
names(umaps_ifn_heterogeneity_list) <- samples
jpeg(paste0(path2figures, "/umap_heterogeneity_only_ifn_genes_negative.jpg"))
gridExtra::grid.arrange(
    grobs=umaps_ifn_heterogeneity_unlisted, 
    ncol=4
)
dev.off()


pdf(paste0(path2figures, "/umap_heterogeneity_only_ifn_genes_negative.pdf"),
    width=18, height = 10)
gridExtra::grid.arrange(
    grobs=umaps_ifn_heterogeneity_unlisted, 
    ncol=4
)
dev.off()




## NFkB signature

## Clustering
seurat_negative <- subset(seurat_2, infected == FALSE)
samples <- unique(seurat_2$orig.ident)
umaps_ifn_heterogeneity_list <- lapply(
    sort(samples), function(samp){
            plot_umap_heterogeneity(seurat_subset = seurat_negative,
                                    sample = samp,
                                    signature = "nfkb_signature",
                                    infection_status = "negative")
    }
)
umaps_ifn_heterogeneity_unlisted <- unlist(umaps_ifn_heterogeneity_list, recursive = FALSE)
names(umaps_ifn_heterogeneity_list) <- samples
jpeg(paste0(path2figures, "/umap_heterogeneity_only_nfkb_genes_negative.jpg"))
gridExtra::grid.arrange(
    grobs=umaps_ifn_heterogeneity_unlisted, 
    ncol=4
)
dev.off()


pdf(paste0(path2figures, "/umap_heterogeneity_only_nfkb_genes_negative.pdf"),
    width=18, height = 10)
gridExtra::grid.arrange(
    grobs=umaps_ifn_heterogeneity_unlisted, 
    ncol=4
)
dev.off()