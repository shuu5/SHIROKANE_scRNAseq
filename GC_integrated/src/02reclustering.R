# libraryの読み込み
library(Seurat)
library(ggplot2)
library(reticulate)

# leidenクラスタリング用に使用するpythonを明示的に指定
use_python("/usr/local/package/python/3.12.0/bin/python", required = TRUE)

# Seuratオブジェクトの読み込み
umap_object <- readRDS("data/processed/umap_object.rds")

# 再クラスタリング
res <- 0.2
umap_object <- FindClusters(umap_object, resolution = res, algorithm = 4, method = "igraph")

# UMAPをプロット
DimPlot(umap_object, reduction = "umap_harmony")
ggsave(file = paste0("output/plot/umap/umap_cluster", res, ".png"), width = 10, height = 10)


# 計算後のseurat_objectを保存
saveRDS(umap_object, file = "data/processed/umap_object.rds")
