# libraryの読み込み
library(Seurat)
library(reticulate)
library(tidyverse)

# leidenクラスタリング用に使用するpythonを明示的に指定
use_python("/usr/local/package/python/3.12.0/bin/python", required = TRUE)

# umap_objectを読み込み
umap_object <- readRDS("data/processed/umap_object.rds")

# celltypeごとのsubclusteringを行う
celltypes <- unique(umap_object$celltype)  # celltypeのリストを取得
subclust_list <- list()  # Seuratオブジェクトを保存するリスト
for (celltype in celltypes) {
    print(paste0("celltype: ", celltype, " is subsetting..."))
    subset_object <- subset(umap_object, idents = celltype)  # celltypeごとにsubsetを作成
    subclust_list[[celltype]] <- subset_object  # リストに保存
}

# 各celltypeごとに再度UMAPを計算してからクラスタリング
for (celltype in names(subclust_list)) {
    print(paste0("celltype: ", celltype, " is running UMAP..."))

    # UMAPを実行
    subclust_list[[celltype]] <- RunUMAP(
        subclust_list[[celltype]], dims = 1:30,
        reduction = "harmony",
        reduction.name = "umap_harmony"
    )
    
    # クラスタリング
    print(paste0("celltype: ", celltype, " is running clustering..."))
    subclust_list[[celltype]] <- FindNeighbors(subclust_list[[celltype]], dims = 1:30)
    subclust_list[[celltype]] <- FindClusters(subclust_list[[celltype]], resolution = 0.5, algorithm = 4, method = "igraph")  
} 

# 各celltypeのオブジェクトに対してUMAPをプロット
for (celltype in names(subclust_list)) {
    print(paste0("celltype: ", celltype, " is plotting UMAP..."))

    DimPlot(subclust_list[[celltype]], reduction = "umap_harmony")
    ggsave(file = paste0("output/plot/umap/subclusters/umap_", celltype, ".png"), width = 10, height = 10)

    DimPlot(subclust_list[[celltype]], reduction = "umap_harmony", group.by = "tissue")
    ggsave(file = paste0("output/plot/umap/subclusters/umap_tissue_", celltype, ".png"), width = 10, height = 10)

    DimPlot(subclust_list[[celltype]], reduction = "umap_harmony", group.by = "batch")
    ggsave(file = paste0("output/plot/umap/subclusters/umap_batch_", celltype, ".png"), width = 10, height = 10)

    DimPlot(subclust_list[[celltype]], reduction = "umap_harmony", group.by = "batch", split.by = "tissue")
    ggsave(file = paste0("output/plot/umap/subclusters/umap_batch_tissue_", celltype, ".png"), width = 15, height = 10)
}

# normalizationに成功したsubclusteringしたオブジェクトをそれぞれrdsに保存
for (celltype in names(subclust_list)) {
    print(paste0("celltype: ", celltype, " is saving subclustering object..."))
    saveRDS(subclust_list[[celltype]], file = paste0("data/processed/subclusters/", celltype, "_subclust.rds"))
}