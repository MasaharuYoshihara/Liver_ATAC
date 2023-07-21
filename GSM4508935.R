# Domcke et al. Science 2020
# GSE149683

library(dplyr)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v75)
library(ggplot2)

seurat_object <- readRDS("GSM4508935_liver_filtered.seurat.RDS")

#head(seurat_object[[]])

DimPlot(seurat_object, reduction = "umap", raster = FALSE, 
        group.by = "cell_type")



png(file = "UMAP_NOTCH_230721.png", width = 520, height = 350)
FeaturePlot(
  object = seurat_object,
  features = c('NOTCH1','NOTCH2','NOTCH3','NOTCH4'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  raster = FALSE
)
dev.off()

png(file = "Dot_NOTCH_230721.png", width = 540, height = 350)
DotPlot(
  object = seurat_object,
  features = c('NOTCH1','NOTCH2','NOTCH3','NOTCH4'),
  group.by = "cell_type"
)
dev.off()



png(file = "UMAP_DSL_230721.png", width = 520, height = 350)
FeaturePlot(
  object = seurat_object,
  features = c('JAG1','JAG2','DLL1','DLL3','DLL4'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  raster = FALSE
)
dev.off()

png(file = "Dot_DSL_230721.png", width = 540, height = 350)
DotPlot(
  object = seurat_object,
  features = c('JAG1','JAG2','DLL1','DLL3','DLL4'),
  group.by = "cell_type"
)
dev.off()


subset_94 <- seurat_object[ , rownames(seurat_object@meta.data)[seurat_object@meta.data$day_of_pregnancy == "94"]]
subset_94

subset_110 <- seurat_object[ , rownames(seurat_object@meta.data)[seurat_object@meta.data$day_of_pregnancy == "110"]]
subset_110

subset_115 <- seurat_object[ , rownames(seurat_object@meta.data)[seurat_object@meta.data$day_of_pregnancy == "115"]]
subset_115

subset_120 <- seurat_object[ , rownames(seurat_object@meta.data)[seurat_object@meta.data$day_of_pregnancy == "120"]]
subset_120

subset_122 <- seurat_object[ , rownames(seurat_object@meta.data)[seurat_object@meta.data$day_of_pregnancy == "122"]]
subset_122

png(file = "dimplot_94_230721.png", width = 550, height = 350)
DimPlot(subset_94, reduction = "umap", raster = FALSE, 
        group.by = "cell_type")
dev.off()

png(file = "dimplot_110_230721.png", width = 550, height = 350)
DimPlot(subset_110, reduction = "umap", raster = FALSE, 
        group.by = "cell_type")
dev.off()

png(file = "dimplot_115_230721.png", width = 550, height = 350)
DimPlot(subset_115, reduction = "umap", raster = FALSE, 
        group.by = "cell_type")
dev.off()

png(file = "dimplot_120_230721.png", width = 550, height = 350)
DimPlot(subset_120, reduction = "umap", raster = FALSE, 
        group.by = "cell_type")
dev.off()

png(file = "dimplot_122_230721.png", width = 550, height = 350)
DimPlot(subset_122, reduction = "umap", raster = FALSE, 
        group.by = "cell_type")
dev.off()

