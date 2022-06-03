library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
set.seed(1234)

args = commandArgs(TRUE)

cur_sample <- args[1]
sample_folder <- args[2]

print(cur_sample)
print(sample_folder)

pbmc.data <- Read10X(data.dir= paste0(sample_folder, "/filtered_feature_bc_matrix/"))

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "infant_flu_vac_RNA", assay = "RNA", min.cells = 3, min.features = 200)

pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

pdf(paste0("VlnPlot_", cur_sample, ".pdf"), width = 15, height = 10)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


pbmc <- subset(pbmc, subset = nFeature_RNA > 100 & nFeature_RNA < 9000)
pbmc

pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# pbmc <- RunUMAP(pbmc, dims = 1:10)

metadata <- readRDS("Infant_Flu_Vac_Metadata.Rds")

pbmc$Sample <- cur_sample

# t.metadata <- pbmc@meta.data

# t.metadata <- merge(t.metadata, metadata, by.x = "Sample", by.y = "DonorID", sort =F)

# pbmc <- AddMetaData(pbmc, t.metadata)

pbmc@meta.data[,colnames(metadata)] <- metadata[cur_sample,,drop=T]

saveRDS(pbmc, file = paste0(cur_sample, ".raw.seurat.Rds"))