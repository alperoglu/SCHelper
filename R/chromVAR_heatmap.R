library(tidyverse)
library(JASPAR2020)
library(TFBSTools)
library(universalmotif)
library(dendextend)
library(grid)
library(gridExtra)
library(ggside)
library(patchwork)

# read in the differential motifs dataframe from Signac Vignette
chromVAR.signature.motifs <- readRDS("~/Desktop/Stitzel_Lab_Work/Agg16_Islets/data/chromVAR.signature.motifs.Rds")

# create a motif database in order to cluster motifs
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# convert motif ids to motif names using database
chromVAR.signature.motifs <- chromVAR.signature.motifs %>%
  rowwise() %>%
  mutate(motif = pfm[[gene]]@name) %>%
  as.data.frame()

# select the top motifs you want to plot for each group
top_motifs <- chromVAR.signature.motifs %>%
  dplyr::group_by(cluster) %>%
  filter(cluster %in% c("Alpha", "Beta", "Delta", "Ductal")) %>%
  slice_min(order_by = p_val_adj, n = 20, with_ties = F) %>%
  as.data.frame() %>%
  dplyr::select(motif,gene)

top_motifs <- top_motifs %>%
  distinct()

# get the top motifs as Universal Motif objects so that we can filter and cluster them
motifs <- universalmotif::convert_motifs(pfm[top_motifs$gene])

# just keep human motifs
motifs <- filter_motifs(motifs, organism = "Homo sapiens")

# plots a hierarchical clustering of the motifs to give an idea
motif_tree(motifs, labels = "name", layout = "rectangular", linecol = "organism")

# find euclidian distance of motifs
comp.matrix <- universalmotif::compare_motifs(motifs, method = "EUCL")

# convert to suitable format for dhc
comp.dist <- as.dist(comp.matrix)

# hierarchical clustering of the motifs for pheatmap
dhc <- hclust(comp.dist)

# just get the top motifs from results and also convert adjusted p value to -log for plotting
df <- chromVAR.signature.motifs %>%
  mutate(cluster = factor(cluster, levels = c("Alpha", "Beta", "Delta", "Gamma", "Ductal", "Acinar", "Endothelial", "Stellate", "Stellate_activated", "Immune"))) %>%
  mutate(Log_Adj_P = -log10(p_val_adj)) %>%
  mutate(Log_Adj_P = ifelse(Log_Adj_P > 400, 400, Log_Adj_P)) %>%
  filter(motif %in% dhc$labels) %>%
  as.data.frame()

# convert data frame to matrix for pheatmap
a <- df[,c("cluster","motif","Log_Adj_P")] %>%
  distinct(cluster, motif, .keep_all = T) %>%
  spread(key = "cluster", value = "Log_Adj_P")

rownames(a) <- a$motif

# pheatmap. a[dhc$labels,-1] is really important as the clustering will not be correct if motifs are not given in the same order as dhc object. Adds stars to boxes based on thresholds
p <- pheatmap::pheatmap(a[dhc$labels,-1], 
                        cluster_cols = F,
                        cluster_rows = as.hclust(dhc),
                        cutree_rows = 24, scale = "none", na_col = "grey",
                        color = colorRampPalette(c("lightblue1","#56B1F7", "#132B43"))(50),
                        angle_col = 45, 
                        display_numbers = matrix(ifelse(a[dhc$labels,-1] > 30, "***", ifelse(a[dhc$labels,-1] > 10, "**",ifelse(a[dhc$labels,-1] > 1.3, "*","" ) )), nrow(a[dhc$labels,-1])),
                        fontsize_number = 12,
                        show_colnames = T, legend = T)


