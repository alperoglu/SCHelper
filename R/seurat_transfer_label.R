library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
set.seed(1234)

args = commandArgs(TRUE)

pbmc_file <- args[1]
ref_file <- args[2]
ref_col_name <- args[3]
out_file <- args[4]

pbmc <- readRDS(pbmc_file)
pbmc_rna <- readRDS(ref_file)

transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  dims = 1:30,
  reference.reduction = "pca"
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna@meta.data[,ref_col_name],
  dims = 1:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

saveRDS(pbmc, file = out_file)