# libraryの読み込み
library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(harmony)
library(reticulate)
library(tidyverse)

# leidenクラスタリング用に使用するpythonを明示的に指定
use_python("/usr/local/package/python/3.12.0/bin/python", required = TRUE)

# グローバル変数の最大サイズを明示的に変更
options(future.globals.maxSize = 2e+09)  # 2GBに設定




# Seuratオブジェクトの読み込み
merged_object <- readRDS("data/processed/merged_seurat.rds")

# # データのサンプル数を1/10にする
# merged_object <- merged_object[, sample(ncol(merged_object), size = floor(ncol(merged_object) / 10))]

#meta.dataの修正
# Laurenの-をNon_Tumorに変更
merged_object@meta.data <- merged_object@meta.data %>%
    mutate(Lauren = ifelse(Lauren == "-", "Non_Tumor", Lauren))

# meta.dataにgroup(Non_Tumor, Tumor)を追加(tissueがNormal or Paratumorの場合はNon_Tumor, それ以外はTumor)
merged_object@meta.data$group <- ifelse(merged_object@meta.data$tissue == "Normal" | merged_object@meta.data$tissue == "Paratumor", "Non_Tumor", "Tumor")




# フィルタリング
merged_object <- subset(
    merged_object, 
    subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10
    )

# batchでlayerを分ける
merged_object[["RNA"]] <- split(merged_object[["RNA"]], f=merged_object@meta.data$batch)

# SCTransform
merged_object <- SCTransform(
    merged_object, 
    min_cells = 3, 
    ncells = ncol(merged_object)
    )

# PCA
merged_object <- RunPCA(merged_object, npcs = 30)

# Harmonyの実行
merged_object <- IntegrateLayers(
    object = merged_object,
    method = HarmonyIntegration,
    normalization.method = "SCT",
    orig.reduction = "pca",
    new.reduction = "harmony"
)

# UMAPを実行
merged_object <- RunUMAP(
    merged_object, dims = 1:30,
    reduction = "harmony",
    reduction.name = "umap_harmony"
    )

# クラスタリング
merged_object <- FindNeighbors(merged_object, dims = 1:30)
merged_object <- FindClusters(merged_object, resolution = 0.1, algorithm = 4, method = "igraph")

# 確認のためにUMAPをプロット
DimPlot(merged_object, reduction = "umap_harmony")
ggsave(file = "output/plot/umap/umap.png", width = 10, height = 10)

DimPlot(merged_object, reduction = "umap_harmony", group.by = "tissue")
ggsave(file = "output/plot/umap/umap_tissue.png", width = 10, height = 10)

DimPlot(merged_object, reduction = "umap_harmony", group.by = "batch")
ggsave(file = "output/plot/umap/umap_batch.png", width = 10, height = 10)

DimPlot(merged_object, reduction = "umap_harmony", group.by = "batch", split.by = "tissue")
ggsave(file = "output/plot/umap/umap_batch_tissue.png", width = 15, height = 10)

# RNA assayのlayerをJoin
DefaultAssay(merged_object) <- "RNA"
merged_object <- JoinLayers(merged_object)
DefaultAssay(merged_object) <- "SCT"

# 計算後のseurat_objectを保存
saveRDS(merged_object, file = "data/processed/umap_object.rds")