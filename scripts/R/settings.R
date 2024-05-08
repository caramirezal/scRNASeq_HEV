## General settings for the scRNASeq analysis

library(RColorBrewer)


##----------------------------------------------------------------------------------
## Definition of the colors

## Samples in the 2nd sequencing
colors_seq2_samples <- c('MockE'='#00cc00',
                           'MockL'='#004d00',
                           'WTE'='#e68a00',
                           'WTEl'='#e60000',
                           'WTL'='#800000',
                           'orf2E1'='#1a8cff',
                           'orf2E2'='#003366')





colors_clusters <- RColorBrewer::brewer.pal(name="Set2", n=8)
cluster_names <- c("C0", "C1", "C2", "C3")
color_clusters <- c("gray50", colors_clusters[1:3])
names(color_clusters) <- cluster_names


colors_infection = c(`TRUE`="red", `FALSE`="gray90")