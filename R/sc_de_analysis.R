library(Seurat)
library(tidyverse)

# pbmc <- readRDS("All.Clusters.seurat.Rds")
pbmc <- readRDS("All.Clusters.v2.seurat.Rds")
DefaultAssay(pbmc)

DefaultAssay(pbmc) <- "RNA"

# pbmc <- subset(pbmc, DonorID %in% c("HPIR-101-6", "HPIR-101-7A"), invert =T)


celltypes <- pbmc$SubAnnotations_v3 %>% unique
Idents(pbmc) <- "SubAnnotations_v3"

vac.markers.list <- lapply(celltypes, function(type){
  vac.markers <- FindMarkers(pbmc, ident.1 = "Post", ident.2 = "Pre", group.by = "Day_2", subset.ident=type, test.use = "negbinom", logfc.threshold = 0, min.pct = 0.1, latent.vars = c("Sex", "Processing.Date", "Race"))
  return(vac.markers)
})

names(vac.markers.list) <- celltypes

saveRDS(vac.markers.list, file = "Vac.Marker.List.new.Annot.Rds")