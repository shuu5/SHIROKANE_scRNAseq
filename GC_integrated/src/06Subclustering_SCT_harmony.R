# libraryの読み込み
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(harmony)
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
    print(paste0("celltype: ", celltype, " is processing..."))
    subset_object <- subset(umap_object, idents = celltype)  # celltypeごとにsubsetを作成
    subclust_list[[celltype]] <- subset_object  # リストに保存
}

# すべてのオブジェクトのデフォルトアッセイをRNAに設定
for (celltype in names(subclust_list)) {
    DefaultAssay(subclust_list[[celltype]]) <- "RNA"  # デフォルトアッセイをRNAに設定
}

# 各celltypeごとにbatchでlayerを分ける
for (celltype in names(subclust_list)) {
    print(paste0("celltype: ", celltype))

    # normalization_successfulをFALSEに設定
    subclust_list[[celltype]]$normalization_successful <- FALSE

    # 細胞が十分なバッチのみ処理
    limit_cells <- 10
    batch_counts <- table(subclust_list[[celltype]]@meta.data$batch)
    valid_batches <- names(batch_counts[batch_counts > limit_cells])
    
    # 2以下の細胞数のバッチを削除
    removed_batches <- names(batch_counts[batch_counts <= limit_cells])
    total_removed_cells <- sum(batch_counts[batch_counts <= limit_cells])  # 削除した細胞の合計数
    if (length(removed_batches) > 0) {
        print(paste0("Removed batches: ", paste(removed_batches, collapse = ", "), 
                     " with a total of ", total_removed_cells, " cells removed."))
    }
    
    # valid_batchesを使ってsubclust_listをフィルタリング
    Idents(subclust_list[[celltype]]) <- subclust_list[[celltype]]@meta.data$batch
    subclust_list[[celltype]] <- subset(subclust_list[[celltype]], 
                                         idents = valid_batches)
    
    # ここでvalid_batchesが空でない場合に処理を続行
    if (length(valid_batches) > 0) {
        # RNA assayのlayerをsplit
        print(paste0("celltype: ", celltype, " is splitting RNA layer..."))
        subclust_list[[celltype]][["RNA"]] <- split(
            subclust_list[[celltype]][["RNA"]], 
            f = subclust_list[[celltype]]@meta.data$batch
        )
        
        # SCTransform
        subclust_list[[celltype]] <- SCTransform(
            subclust_list[[celltype]], 
            min_cells = 3, 
            ncells = ncol(subclust_list[[celltype]])
        )
        
        # PCA
        subclust_list[[celltype]] <- RunPCA(subclust_list[[celltype]], npcs = 30)
        
        # Harmonyの実行
        subclust_list[[celltype]] <- IntegrateLayers(
            object = subclust_list[[celltype]],
            method = HarmonyIntegration,
            normalization.method = "SCT",
            orig.reduction = "pca",
            new.reduction = "harmony"
        )
        
        # UMAPを実行
        subclust_list[[celltype]] <- RunUMAP(
            subclust_list[[celltype]], dims = 1:30,
            reduction = "harmony",
            reduction.name = "umap_harmony"
        )
        
        # クラスタリング
        subclust_list[[celltype]] <- FindNeighbors(subclust_list[[celltype]], dims = 1:30)
        subclust_list[[celltype]] <- FindClusters(subclust_list[[celltype]], resolution = 0.5, algorithm = 4, method = "igraph")
        
        # 正規化処理が完了した場合はフラグを更新
        subclust_list[[celltype]]$normalization_successful <- TRUE
    } else {
        warning(paste0("celltype: ", celltype, " has no valid batches and will be skipped."))
    }
    
    # RNA assayのlayerをJoin
    DefaultAssay(subclust_list[[celltype]]) <- "RNA"
    subclust_list[[celltype]] <- JoinLayers(subclust_list[[celltype]])
    DefaultAssay(subclust_list[[celltype]]) <- "SCT"
}


# normalizationに成功したcelltypeのオブジェクトに対してUMAPをプロット
for (celltype in names(subclust_list)) {
    if (subclust_list[[celltype]]$normalization_successful[1]) {
        print(paste0("celltype: ", celltype, " is plotting UMAP 1/4..."))
        DimPlot(subclust_list[[celltype]], reduction = "umap_harmony")
        ggsave(file = paste0("output/plot/umap/subclusters/umap_", celltype, ".png"), width = 10, height = 10)

        print(paste0("celltype: ", celltype, " is plotting UMAP 2/4..."))
        DimPlot(subclust_list[[celltype]], reduction = "umap_harmony", group.by = "tissue")
        ggsave(file = paste0("output/plot/umap/subclusters/umap_tissue_", celltype, ".png"), width = 10, height = 10)

        print(paste0("celltype: ", celltype, " is plotting UMAP 3/4..."))
        DimPlot(subclust_list[[celltype]], reduction = "umap_harmony", group.by = "batch")
        ggsave(file = paste0("output/plot/umap/subclusters/umap_batch_", celltype, ".png"), width = 10, height = 10)

        print(paste0("celltype: ", celltype, " is plotting UMAP 4/4..."))
        DimPlot(subclust_list[[celltype]], reduction = "umap_harmony", group.by = "batch", split.by = "tissue")
        ggsave(file = paste0("output/plot/umap/subclusters/umap_batch_tissue_", celltype, ".png"), width = 15, height = 10)
    }
}

# normalizationに成功したsubclusteringしたオブジェクトをそれぞれrdsに保存
for (celltype in names(subclust_list)) {
    if (subclust_list[[celltype]]$normalization_successful[1]) {
        print(paste0("celltype: ", celltype, " is saving..."))
        saveRDS(subclust_list[[celltype]], file = paste0("data/processed/subclusters/", celltype, "_subclust.rds"))
    }
}