# libraryの読み込み
library(Seurat)
library(reticulate)
library(tidyverse)

recluster_cells <- function(cell_type, res = 0.5, algorithm = 4, method = "igraph") {
    print(paste0("celltype: ", cell_type, " is reclustering..."))

    # RDSファイルのパスを作成
    file_path <- paste0("data/processed/subclusters/", cell_type, "_subclust.rds")
    
    # Seuratオブジェクトを読み込み
    seurat_obj <- readRDS(file_path)
    
    # 再クラスタリング
    seurat_obj <- FindClusters(seurat_obj, resolution = res, algorithm = algorithm, method = method)

    # UMAP plot
    print(paste0("celltype: ", cell_type, " is plotting UMAP..."))

    DimPlot(seurat_obj, reduction = "umap_harmony")
    ggsave(file = paste0("output/plot/umap/subclusters/umap_", cell_type, ".png"), width = 10, height = 10)

    DimPlot(seurat_obj, reduction = "umap_harmony", group.by = "tissue")
    ggsave(file = paste0("output/plot/umap/subclusters/umap_tissue_", cell_type, ".png"), width = 10, height = 10)

    DimPlot(seurat_obj, reduction = "umap_harmony", group.by = "batch")
    ggsave(file = paste0("output/plot/umap/subclusters/umap_batch_", cell_type, ".png"), width = 10, height = 10)

    DimPlot(seurat_obj, reduction = "umap_harmony", group.by = "batch", split.by = "tissue")
    ggsave(file = paste0("output/plot/umap/subclusters/umap_batch_tissue_", cell_type, ".png"), width = 15, height = 10)
    
    # 結果を保存
    saveRDS(seurat_obj, file = file_path)
    
    return(seurat_obj)
}



# data/processed/subclusters/ディレクトリのファイル名を表示
list.files("data/processed/subclusters/", pattern = "*.rds", full.names = TRUE)

# 再クラスタリングを実行
recluster_cells("Epithelial_cells", res = 0.15, algorithm = 4, method = "igraph")
