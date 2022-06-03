library(Seurat)
library(tidyverse)
library(DoubletFinder)

args = commandArgs(TRUE)

cur_sample <- args[1]
sample_folder <- args[2]

print(cur_sample)
print(sample_folder)

# pbmc.data <- Read10X_h5(filename = paste0(sample_folder, "/filtered_feature_bc_matrix.h5"))
pbmc.data <- Read10X(data.dir= paste0(sample_folder, "/filtered_feature_bc_matrix/"))

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
nExp_poi <- round(0.075*nrow(pbmc@meta.data))
pbmc <- doubletFinder_v3(pbmc, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

doublets <- pbmc@meta.data[,5,drop = F]
colnames(doublets) <- "DF_results"

print(table(doublets$DF_results))

write.table(row.names(doublets[doublets$DF_results == "Doublet",,drop = F]),
            quote = F, row.names = F, col.names = F, file = paste0("doublets_df_", cur_sample, ".txt"))

rm(pbmc, pbmc.data)
gc()