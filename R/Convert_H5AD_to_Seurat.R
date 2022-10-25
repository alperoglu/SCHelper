library(anndata)
library(tidyverse)
library(Seurat)
library(Matrix)

args = commandArgs(TRUE)

input_adata <- args[1]
output_seurat <- args[2]

pbmc_adata <- read_h5ad(input_adata)

a <- pbmc_adata$raw$to_adata()

data <- a$to_df() %>% as.matrix()

print(dim(data))

counts = expm1(data)

metadata <- pbmc_adata$obs

print(dim(metadata))

rm(data, a)

total_counts <- pbmc_adata$obs['n_counts'] %>% as.data.frame()

t.counts <- counts * total_counts[,1]

rm(counts, pbmc_adata)

t.counts <- t.counts / 10000

t.counts <- floor(t.counts)

print(dim(t.counts))

pbmc_seurat<- CreateSeuratObject(counts = Matrix(t(t.counts), sparse =T), min.cells = 0, min.features = 0)

pbmc_seurat <- AddMetaData(pbmc_seurat, metadata = metadata)

print(pbmc_seurat)

# data  <- pbmc$obsm['X_pca'][[1]]

# rownames(data) <- pbmc$obs_names
# colnames(data) <- paste0("PC_", 1:ncol(data))
# pbmc_seurat[["pca"]] <- CreateDimReducObject(embeddings = data, assay = "RNA", key = "PC")

# data  <- pbmc_adata$obsm['X_umap'][[1]]

# rownames(data) <- pbmc_adata$obs_names
# colnames(data) <- paste0("UMAP_", 1:ncol(data))

# pbmc_seurat[["umap"]] <- CreateDimReducObject(embeddings = data, assay = "RNA", key = "UMAP")

rm(pbmc_adata, a, data, counts, t.counts, total_counts)

# print("before SCTransform")

# pbmc_seurat <- SCTransform(pbmc_seurat, verbose = T, conserve.memory = T)

saveRDS(pbmc_seurat, file  = output_seurat)
