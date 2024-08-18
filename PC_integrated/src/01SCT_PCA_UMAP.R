# Seuratパッケージの読み込み
library(Seurat)
library(harmony)
library(tidyverse)
library(sctransform)


# フィルタリング
seurat_obj_pdac <- subset(seurat_obj_pdac, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# batchでlayerを分ける
seurat_obj_pdac[["RNA"]] <- split(seurat_obj_pdac[["RNA"]], f = seurat_obj_pdac$batch)

# SCTransformを使用した正規化と変動遺伝子の特定
seurat_obj_pdac <- SCTransform(seurat_obj_pdac)

# PCAの実行
seurat_obj_pdac <- RunPCA(seurat_obj_pdac, npcs = 30)

# Harmonyの実行
seurat_obj_pdac <- IntegrateLayers(
    object = seurat_obj_pdac,
    method = HarmonyIntegration,
    normalization.method = "SCT",
    orig.reduction = "pca",
    new.reduction = "harmony"
)

# UMAPの実行
seurat_obj_pdac <- RunUMAP(
    seurat_obj_pdac, dims = 1:30,
    reduction = "harmony",
    reduction.name = "umap_harmony"
    )


# クラスタリング
seurat_obj_pdac <- FindNeighbors(seurat_obj_pdac, dims = 1:30)
seurat_obj_pdac <- FindClusters(
    seurat_obj_pdac, 
    resolution = 0.1,
    algorithm = 4,
    method = "igraph"
    )

# UMAPによる可視化
DimPlot(seurat_obj_pdac, reduction = "umap_harmony", group.by = "batch")
ggsave("output/plot/UMAP_batch.png")

# 結果の表示
print(seurat_obj_pdac)

# seurat_obj_pdacを保存
saveRDS(seurat_obj_pdac, "data/processed/seurat_obj_pdac.rds")
