# libraryの読み込み
library(Seurat)
library(magrittr)
library(glmGamPoi)
library(patchwork)
library(SingleR)
library(SingleCellExperiment)
library(celldex)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)

# folder name
folder_name <- "HB_GSE186975"

# seurat_objectを読み込み
seurat_object <- readRDS(paste0("data/", folder_name, "/merged_seurat_object.rds"))

### 統合後のデータのPCA、クラスタリング、可視化
# デフォルトアッセイの指定
DefaultAssay(seurat_object) <- "RNA"

# 変数特徴量の特定
seurat_object <- FindVariableFeatures(seurat_object)

# PCAの実行
seurat_object <- RunPCA(seurat_object, npcs = 30, verbose = FALSE)

# 寄与率の確認
ElbowPlot(seurat_object, ndims = 30)

# FindNeighbors
seurat_object <- FindNeighbors(seurat_object, dims = 1:30)

# FindClusters
seurat_object <- FindClusters(seurat_object, resolution = 0.5)

# UMAPを実行
seurat_object <- RunUMAP(
  seurat_object,
  reduction = "pca",
  dims = 1:30,
  umap.method = "umap-learn",
  metric = "correlation"
)

# UMAPをプロット
DimPlot(seurat_object, reduction = "umap", group.by = "batch", split.by = "tissue", label = TRUE, pt.size = 0.5)


# 計算後のseurat_objectを保存
saveRDS(seurat_object, paste0("data/", folder_name, "/seurat_object_UMAP.rds"))
