## Heatmap of gene expression

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


path2figures <- paste0(path2project, '/figures/03_gex_heatmaps')
if ( ! dir.exists(path2figures)) {
    dir.create(path2figures)
}

path2analysis <- paste0(path2project, '/analysis/03_gex_heatmaps')
if ( ! dir.exists(path2analysis)) {
    dir.create(path2analysis)
}


## Reading seurat objects
seurat_2 <- readRDS('/media/ag-cherrmann/projects/24_scRNA_Loan/analysis/seuratIntegrated_runsJan2024_filtprocd.rds')
seurat_2



## Loading signatures
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




## Annotating seurat object with clusters
path2clusterAnns <- paste0(path2project, "/analysis/02_bystander_cells/")
anns_ifn_negative <- read_tsv(paste0(path2clusterAnns, 
                                     "/mapping_renamed_clusters_negative_ifn.tsv.gz"))
anns_ifn_all <- read_tsv(paste0(path2clusterAnns, 
                                "/mapping_renamed_clusters_all_ifn.tsv.gz"))                      
anns_nfkb_negative <- read_tsv(paste0(path2clusterAnns, 
                                     "/mapping_renamed_clusters_negative_nfkb.tsv.gz"))
anns_nfkb_all <- read_tsv(paste0(path2clusterAnns, 
                                "/mapping_renamed_clusters_all_nfkb.tsv.gz"))   


## Adding ifn and nfkb siganture expression
path2auc <- paste0(path2project, "/analysis/01_ifn_and_nfkb_responses/")
signatures_auc_df <- read_tsv(paste0(path2auc, "/ifn_nfkb_aucell.tsv.gz"))
all(signatures_auc_df$barcode==colnames(seurat_2))
seurat_2$"ifn_signature" <- signatures_auc_df$"ifn_signature"
seurat_2$"nfkb_signature" <- signatures_auc_df$"nfkb_signature"



sample_name <- "WTL"
#signature_name <- "ifn_signature"
infection_status <- "all"


    signature_df = filter(signature_df, ORF2_log < quantile(ORF2_log, probs = 0.95, na.rm = TRUE))

plot_heatmap <- function(
    sample_name,
    signature_name
){
    signature_df = FetchData(subset(seurat_2, 
                                    orig.ident==sample_name &
                                      ifn_signature < quantile(
                                                ifn_signature, 
                                                probs = 0.95,
                                                 na.rm = TRUE) &
                                      nfkb_signature < quantile(
                                                nfkb_signature, 
                                                probs = 0.95,
                                                 na.rm = TRUE) &
                                      ORF2 < quantile(ORF2, probs = 0.95)
                                     ), 
                        vars=c(signatures_list[signature_name][[1]], 
                                "ORF2", "ifn_signature",
                                "nfkb_signature", "infected"), 
                        layer="data")
    signature_df = mutate(signature_df, ORF2_log=log10(ORF2+1))
    signature_df$"cluster_renamed" <- plyr::mapvalues(
    x = rownames(signature_df),
    from =  anns_ifn_all$"barcode",
    to = anns_ifn_all$"cluster_renamed"
    )   
    #table(signature_df$cluster_renamed)
    mtx = select(signature_df, 
                -ORF2, 
                -ORF2_log, 
                -ifn_signature, 
                -cluster_renamed, 
                -nfkb_signature,
                -infected) %>%
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
                    cluster_renamed,
                    ifn_signature,
                    nfkb_signature,
                    infected)
    min <- min(pull(mdata, ifn_signature))
    max = max(pull(mdata, ifn_signature))
    mean_ifn <- mean(pull(mdata, ifn_signature))
    col_fun_ifn = circlize::colorRamp2(breaks=c(min,
                                                mean_ifn,
                                                max), 
                                        c("blue",
                                        "yellow",
                                        "red"))
    min <- min(pull(mdata, nfkb_signature))
    max = max(pull(mdata, nfkb_signature))
    mean_nfkb <- mean(pull(mdata, nfkb_signature))
    col_fun_nfkb = circlize::colorRamp2(breaks=c(min,
                                                mean_nfkb,
                                                max), 
                                        c("blue",
                                        "yellow",
                                        "red"))
    min <- min(pull(mdata, ORF2_log))
    max = max(pull(mdata, ORF2_log))
    mean_orf2 <- mean(pull(mdata, ORF2_log))
    col_fun_orf2 = circlize::colorRamp2(breaks=c(min,
                                                mean_orf2,
                                                max), 
                                        c("blue",
                                        "yellow",
                                        "red"))
    col_ann <- HeatmapAnnotation(ifn_signature=mdata$"ifn_signature",
                                nfkb_signature=mdata$"nfkb_signature",
                                cluster_renamed=mdata$"cluster_renamed",
                                ORF2=mdata$"ORF2_log",
                                col = list(ifn_signature=col_fun_ifn,
                                            nfkb_signature=col_fun_nfkb,
                                            cluster_renamed=color_clusters,
                                            ORF2=col_fun_orf2
                                            ))

    Heatmap(mtx,
            top_annotation = col_ann,
            show_column_names = FALSE,
            column_split = mdata$infected,
            column_title = sample_name)

}



## Stupid vscode compiler is not working properly
## I need to split the function in two
plot_heatmap_mock <- function(
    sample_name,
    signature_name
){
    signature_df = FetchData(subset(seurat_2, 
                                    orig.ident==sample_name &
                                      ifn_signature < quantile(
                                                ifn_signature, 
                                                probs = 0.95,
                                                 na.rm = TRUE) &
                                      nfkb_signature < quantile(
                                                nfkb_signature, 
                                                probs = 0.95,
                                                 na.rm = TRUE) 
                                     ), 
                        vars=c(signatures_list[signature_name][[1]], 
                                "ORF2", "ifn_signature",
                                "nfkb_signature"), 
                        layer="data")
    signature_df$"cluster_renamed" <- plyr::mapvalues(
    x = rownames(signature_df),
    from =  anns_ifn_all$"barcode",
    to = anns_ifn_all$"cluster_renamed"
    )   
    #table(signature_df$cluster_renamed)
    mtx = select(signature_df, 
                -ORF2, 
                -ifn_signature, 
                -cluster_renamed, 
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
                    cluster_renamed,
                    ifn_signature,
                    nfkb_signature)
    min <- min(pull(mdata, ifn_signature))
    max = max(pull(mdata, ifn_signature))
    mean_ifn <- mean(pull(mdata, ifn_signature))
    col_fun_ifn = circlize::colorRamp2(breaks=c(min,
                                                mean_ifn,
                                                max), 
                                        c("blue",
                                        "yellow",
                                        "red"))
    min <- min(pull(mdata, nfkb_signature))
    max = max(pull(mdata, nfkb_signature))
    mean_nfkb <- mean(pull(mdata, nfkb_signature))
    col_fun_nfkb = circlize::colorRamp2(breaks=c(min,
                                                mean_nfkb,
                                                max), 
                                        c("blue",
                                        "yellow",
                                        "red"))
    col_ann <- HeatmapAnnotation(ifn_signature=mdata$"ifn_signature",
                                nfkb_signature=mdata$"nfkb_signature",
                                cluster_renamed=mdata$"cluster_renamed",
                                col = list(ifn_signature=col_fun_ifn,
                                            nfkb_signature=col_fun_nfkb,
                                            cluster_renamed=color_clusters
                                            ))

    Heatmap(mtx,
            top_annotation = col_ann,
            show_column_names = FALSE,
            column_title = sample_name)

}



samples <- unique(seurat_2$orig.ident) %>% sort()
names(samples) <- samples
heatmaps_ifn_list <- lapply(
    samples, function(samp){
        if ( grepl("Mock", samp )) {
            hm <- plot_heatmap_mock(sample_name = samp,
                         signature_name = "ifn_signature")
        } else {
            hm <- plot_heatmap(sample_name = samp,
                         signature_name = "ifn_signature")
        }
        return(hm)
    }
)




pdf(paste0(path2figures, "/heatmap_gex_ifn.pdf"),
    width = 60, height=10)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 7)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(heatmaps_ifn_list$"MockE", newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(heatmaps_ifn_list$"MockL", newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
draw(heatmaps_ifn_list$"WTE", newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
draw(heatmaps_ifn_list$"WTEl", newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 5))
draw(heatmaps_ifn_list$"WTL", newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 6))
draw(heatmaps_ifn_list$"orf2E1", newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 7))
draw(heatmaps_ifn_list$"orf2E2", newpage = FALSE)
upViewport()

dev.off()



##---------------------------------------------------------------------------
## NFkB 

heatmaps_nfkb_list <- lapply(
    samples, function(samp){
        if ( grepl("Mock", samp )) {
            hm <- plot_heatmap_mock(sample_name = samp,
                         signature_name = "nfkb_signature")
        } else {
            hm <- plot_heatmap(sample_name = samp,
                         signature_name = "nfkb_signature")
        }
        return(hm)
    }
)




pdf(paste0(path2figures, "/heatmap_gex_nfkb.pdf"),
    width = 60, height=10)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 7)))

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(heatmaps_nfkb_list$"MockE", newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(heatmaps_nfkb_list$"MockL", newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
draw(heatmaps_nfkb_list$"WTE", newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
draw(heatmaps_nfkb_list$"WTEl", newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 5))
draw(heatmaps_nfkb_list$"WTL", newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 6))
draw(heatmaps_nfkb_list$"orf2E1", newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 7))
draw(heatmaps_nfkb_list$"orf2E2", newpage = FALSE)
upViewport()

dev.off()
