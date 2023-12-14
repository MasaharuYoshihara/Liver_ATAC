library(Seurat)

rna_D1 <- readRDS("GSE171993_D1.rds")
# 24043 features across 9623 samples
rna_D1_1 <- subset(rna_D1, subset = lib_ID == 1)
# 24043 features across 2405 samples
rna_D1_2 <- subset(rna_D1, subset = lib_ID == 2)
# 24043 features across 7218 samples

DimPlot(rna_D1, reduction = "umap")
DimPlot(rna_D1, reduction = "umap", split.by = "lib_ID")

FeaturePlot(rna_D1, "Pecam1") # endo
FeaturePlot(rna_D1, "Acta2")  # SM
FeaturePlot(rna_D1, "Alb")    # hepatocyte
FeaturePlot(rna_D1, "Epcam")  # cholangiocyte
FeaturePlot(rna_D1, "Gypa")   # erythroid
FeaturePlot(rna_D1, "Lyz2")   # neutrophil 
FeaturePlot(rna_D1, "Nkg7")   # NK/T 
FeaturePlot(rna_D1, "Cd79a")  # B cell

new.cluster.ids <- c("B cell", "Neutro_1", "Ery_1", 
                     "Ery_2", "Neutro_2", "Neutro_3",
                     "Hepatocyte", "Neutro_4", "endo/SM",
                     "Ery_3", "NK/T")
names(new.cluster.ids) <- levels(rna_D1)
rna_D1 <- RenameIdents(rna_D1, new.cluster.ids)
DimPlot(rna_D1, reduction = "umap", label = TRUE)
DimPlot(rna_D1, reduction = "umap", label = TRUE, split.by = "lib_ID")

markers <- c("Pecam1", "Acta2", "Alb", "Epcam", "Gypa", "Lyz2", "Nkg7", "Cd79a")
DotPlot(rna_D1, features = markers)

FeaturePlot(rna_D1, "Jag1")
FeaturePlot(rna_D1, features =  c("Jag2", "Dll1", "Dll3", "Dll4"))
FeaturePlot(rna_D1, features = c("Notch1", "Notch2", "Notch3", "Notch4"))
FeaturePlot(rna_D1, features = c("Hes1", "Hes5", "Hey1", "Hey2"))

D1_Pecam1 <- subset(rna_D1, subset = Pecam1 >0)        # 784 cells
D1_Acta2 <- subset(rna_D1, subset = Acta2 >0)          # 272 cells 
D1_Jag1 <- subset(rna_D1, subset = Jag1 >0)            # 92 cells
D1_Pecam1_Jag1 <- subset(D1_Jag1, subset = Pecam1 >0)  # 24 cells
D1_Acta2_Jag1 <- subset(D1_Jag1, subset = Acta2 >0)    # 10 cells
D1_Pecam1_Acta2 <- subset(rna_D1, subset = Pecam1 > 0 & Acta2 >0) # 22 cells

# count <- Read10X_h5("GSE171993_filtered_feature_bc_matrix.h5")
# count <- CreateSeuratObject(counts = count, min.cells = 3, min.features = 200)
# 
# metadata <- read.csv("GSE171993_metadata.csv")
# rownames(metadata) <- metadata$X
# count <- AddMetaData(object = count, metadata = metadata)
# 
# count[["percent.mt"]] <- PercentageFeatureSet(count, pattern = "^mt-")
# VlnPlot(count, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# count <- subset(count, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25)
# 
# count_D1 <- subset(count, subset = timepoint == "D1")
# 
# count_D1 <- NormalizeData(count_D1, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# count_D1 <- FindVariableFeatures(count_D1, selection.method = "vst", nfeatures = 2000)
# 
# all.genes <- rownames(count_D1)
# count_D1 <- ScaleData(count_D1, features = all.genes)
# 
# count_D1 <- RunPCA(count_D1, features = VariableFeatures(object = count_D1))
# 
# ElbowPlot(count_D1)
# 
# count_D1 <- FindNeighbors(count_D1, dims = 1:15)
# count_D1_res0.2 <- FindClusters(count_D1, resolution = 0.2)
# 
# count_D1_res0.2 <- RunUMAP(count_D1_res0.2, dims = 1:15)
# 
# DimPlot(count_D1_res0.2, reduction = "umap")
# 
# saveRDS(count_D1_res0.2, file = "GSE171993_D1.rds")
# 
# sessionInfo()
# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.2
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Seurat_5.0.1       SeuratObject_5.0.1 sp_2.1-2          
# 
# loaded via a namespace (and not attached):
#   [1] deldir_2.0-2           pbapply_1.7-2         
# [3] gridExtra_2.3          rlang_1.1.2           
# [5] magrittr_2.0.3         RcppAnnoy_0.0.21      
# [7] matrixStats_1.2.0      ggridges_0.5.4        
# [9] compiler_4.3.1         spatstat.geom_3.2-7   
# [11] png_0.1-8              vctrs_0.6.5           
# [13] reshape2_1.4.4         hdf5r_1.3.8           
# [15] stringr_1.5.1          pkgconfig_2.0.3       
# [17] fastmap_1.1.1          ellipsis_0.3.2        
# [19] labeling_0.4.3         utf8_1.2.4            
# [21] promises_1.2.1         ggbeeswarm_0.7.2      
# [23] bit_4.0.5              purrr_1.0.2           
# [25] jsonlite_1.8.8         goftest_1.2-3         
# [27] later_1.3.2            spatstat.utils_3.0-4  
# [29] irlba_2.3.5.1          parallel_4.3.1        
# [31] cluster_2.1.6          R6_2.5.1              
# [33] ica_1.0-3              stringi_1.8.3         
# [35] RColorBrewer_1.1-3     spatstat.data_3.0-3   
# [37] reticulate_1.34.0      parallelly_1.36.0     
# [39] lmtest_0.9-40          scattermore_1.2       
# [41] Rcpp_1.0.11            tensor_1.5            
# [43] future.apply_1.11.0    zoo_1.8-12            
# [45] sctransform_0.4.1      httpuv_1.6.13         
# [47] Matrix_1.6-4           splines_4.3.1         
# [49] igraph_1.6.0           tidyselect_1.2.0      
# [51] rstudioapi_0.15.0      abind_1.4-5           
# [53] spatstat.random_3.2-2  codetools_0.2-19      
# [55] miniUI_0.1.1.1         spatstat.explore_3.2-5
# [57] listenv_0.9.0          lattice_0.22-5        
# [59] tibble_3.2.1           plyr_1.8.9            
# [61] withr_2.5.2            shiny_1.8.0           
# [63] ROCR_1.0-11            ggrastr_1.0.2         
# [65] Rtsne_0.17             future_1.33.0         
# [67] fastDummies_1.7.3      survival_3.5-7        
# [69] polyclip_1.10-6        fitdistrplus_1.1-11   
# [71] pillar_1.9.0           KernSmooth_2.23-22    
# [73] plotly_4.10.3          generics_0.1.3        
# [75] RcppHNSW_0.5.0         ggplot2_3.4.4         
# [77] munsell_0.5.0          scales_1.3.0          
# [79] globals_0.16.2         xtable_1.8-4          
# [81] glue_1.6.2             lazyeval_0.2.2        
# [83] tools_4.3.1            data.table_1.14.10    
# [85] RSpectra_0.16-1        RANN_2.6.1            
# [87] leiden_0.4.3.1         dotCall64_1.1-1       
# [89] cowplot_1.1.1          grid_4.3.1            
# [91] tidyr_1.3.0            colorspace_2.1-0      
# [93] nlme_3.1-164           patchwork_1.1.3       
# [95] beeswarm_0.4.0         vipor_0.4.5           
# [97] cli_3.6.2              spatstat.sparse_3.0-3 
# [99] spam_2.10-0            fansi_1.0.6           
# [101] viridisLite_0.4.2      dplyr_1.1.4           
# [103] uwot_0.1.16            gtable_0.3.4          
# [105] digest_0.6.33          progressr_0.14.0      
# [107] ggrepel_0.9.4          farver_2.1.1          
# [109] htmlwidgets_1.6.4      htmltools_0.5.7       
# [111] lifecycle_1.0.4        httr_1.4.7            
# [113] mime_0.12              bit64_4.0.5           
# [115] MASS_7.3-60  