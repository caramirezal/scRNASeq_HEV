## Analysis of the microarray data in order to infer temporal pseudotime

## 
library(dplyr)
library(ggplot2)
library(factoextra)
library(viridis)
library(ggrepel)



set.seed(333)

rerun <- FALSE

## Working directory
path2project <- '/media/ag-cherrmann/cramirez/scRNASeq_HEV'
setwd(path2project)


## Loading general settings
source(paste0(path2project, '/scripts/R/settings.R'))


path2figures <- paste0(path2project, '/figures/05_bulk_pseudotime')
if ( ! dir.exists(path2figures)) {
    dir.create(path2figures)
}

path2analysis <- paste0(path2project, '/analysis/05_bulk_pseudotime')
if ( ! dir.exists(path2analysis)) {
    dir.create(path2analysis)
}



##----------------------------------------------------------------------------------------
## Ana-Luisa script chunk 


invisible(setwd("/media/ag-cherrmann/lcosta/Interferon/src/microarray_analysis"))
normalised <- read.table("Normalised.txt",header=T)

# filter out non-significant probes
normalised <- subset(normalised, 
                     ((apply(normalised[,c(5,seq(69, 109, 8))] < 0.01, 1, sum)) > 2) &
                       ((apply(normalised[,c(9,seq(73, 113, 8))] < 0.01, 1, sum)) > 2))
normalised_complete <- normalised

normalised <- normalised[,grep("0.noRNA|epo|Symbol|GeneName|EntrezID|Probe", 
                               colnames(normalised))]

# making a df with only intensities (averages)
dmMean <- data.matrix(normalised[,grep(".mean$", colnames(normalised))])




# Save all the names in separate df
# Get just the ones we want (they have dots in them)
names_nr = data.frame(full = colnames(normalised)[grep("\\.", 
                                                       colnames(normalised))])

# Separate into columns
names_wide = tidyr::separate(names_nr, 
                             col = "full", 
                             sep = "\\.",
                             into = c("replicate", "type", "time", "info", "metricx"))
# Try to make the "key"
key = names_wide %>%
  select(replicate, type, time, info) %>%
  distinct()





pca_mod <- prcomp(t(dmMean))  # compute principal components
df_pc <- data.frame(pca_mod$x, 
                    time=key$time , 
                    Genotype = key$type , 
                    replicate = as.factor(key$replicate), 
                    info = key$info  )  # dataframe of principal components

df_pc$time <- factor(df_pc$time,
                     levels = c(0,2,4,6,8,16,24))


##-----------------------------------------------------------------------------------




##-----------------------------------------------------------------------------------
## Adding the pseudotime

ang <- function(x,y) { 
  z <- x + 1i * y
  res <- -140 - Arg(z) / pi * 180
  res %% 360
}



df_pc$"pseudotime" <-  atan2(df_pc$PC1, df_pc$PC2) 
df_pc$"pseudotime" <- abs(min(df_pc$"pseudotime")) + df_pc$"pseudotime" 


## Reproducing Ana's plot
jpeg(paste0(path2figures, "/pca_bulk_microarray.jpg"))
df_pc %>%
    ggplot(aes(x=PC1, y=PC2,
          colour=pseudotime,
          label=time)) +
            geom_text_repel(colour="black", force=4) +
           geom_point(size=5) +
           scale_color_viridis() +
           theme_classic()  +
           xlim(-70, 70) +
           ylim(-35, 30)
dev.off()




##----------------------------------------------------------------------------
## Proyecting samples into clock 

# Sampling 4 points
sample_names <- colnames(dmMean)
sample_indexes <- sample(1:ncol(dmMean), 4, replace = FALSE)
names(sample_indexes) <- sample_names[sample_indexes]
pca_mod <- prcomp(t(dmMean[,-sample_indexes]))  # compute principal components
df_pc <- data.frame(pca_mod$x, 
                    time=key$time[-sample_indexes], 
                    Genotype = key$type[-sample_indexes], 
                    replicate = as.factor(key$replicate[-sample_indexes]), 
                    info = key$info[-sample_indexes] )  # dataframe of principal components

df_pc$time <- factor(df_pc$time,
                     levels = c(0,2,4,6,8,16,24))




## Predicting pseudotime
preds <- predict(pca_mod, newdata=t(dmMean[, sample_indexes]))
preds_df <- as.data.frame(preds)
preds_df <- data.frame(preds_df, 
                    time=key$time[sample_indexes], 
                    Genotype = key$type[sample_indexes], 
                    replicate = as.factor(key$replicate[sample_indexes]), 
                    info = key$info[sample_indexes] )  # dataframe of principal components


df <- rbind(df_pc, preds_df)
df$time <- factor(df$time, levels = c(0,2,4,6,8,16,24))
df$"pseudotime" <-  atan2(df$PC1, df$PC2) 
df$"pseudotime" <- abs(min(df$"pseudotime")) + df$"pseudotime" 



## Prediction subsample
pdf(paste0(path2figures, "/pca_bulk_microarray_predictions.pdf"),
    width = 6, height = 4)
df[!rownames(df) %in% rownames(preds_df), ] %>%
    ggplot(aes(x=PC1, y=PC2,
               colour=pseudotime,
               label=time)) +
           geom_point(size=6, shape=21) +
           geom_point(data=df[rownames(df) %in% rownames(preds_df),],
                        aes(x=PC1, y=PC2,
                            colour=pseudotime,
                            label=time),
                            shape=23,
                            size=6) +
           scale_color_viridis() +
           geom_text_repel(colour="black", force=8) +
           theme_classic()  +
           xlim(-60, 60) +
           ylim(-35, 30)
dev.off()

