library(Seurat)
library(patchwork)

seurat.obj.list <- Sys.glob("*.raw.seurat.Rds")

names(seurat.obj.list) <- gsub(pattern=".raw.seurat.Rds", replacement = "", x = seurat.obj.list)

seurat.obj.list <- lapply(names(seurat.obj.list), function(x){
  obj <- readRDS(seurat.obj.list[x])
  return(obj)
})

library(tidyverse)
names(seurat.obj.list) <- Sys.glob("*.raw.seurat.Rds") %>% gsub(pattern=".raw.seurat.Rds", replacement = "", x = .) %>%gsub(pattern="-", replacement="_", x = .)

features <- SelectIntegrationFeatures(object.list =seurat.obj.list)
pbmc.anchors <-  FindIntegrationAnchors(object.list =seurat.obj.list, anchor.features=features)

saveRDS(pbmc.anchors, file = "Combined.Integration.Anchors.Rds")

pbmc.combined <- IntegrateData(anchorset = pbmc.anchors)

DefaultAssay(pbmc.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
pbmc.combined <- ScaleData(pbmc.combined, verbose = FALSE)
pbmc.combined <- RunPCA(pbmc.combined, npcs = 30, verbose = FALSE)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "pca", dims = 1:30)

saveRDS(pbmc.combined, file = "Combined.raw.seurat.v2.Rds")