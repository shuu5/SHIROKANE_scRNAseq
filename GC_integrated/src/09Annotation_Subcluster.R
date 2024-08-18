library(Seurat)
library(tidyverse)

# data/processed/subclusters/内のsubclusteringしたオブジェクトを読み込み
subclust_list <- lapply(list.files("data/processed/subclusters/", pattern = "*.rds", full.names = TRUE), readRDS)
names(subclust_list) <- gsub(
    "data/processed/subclusters/(.*)_subclust.rds", "\\1", 
    list.files("data/processed/subclusters/", pattern = "*.rds")
)
names(subclust_list) <- gsub("_subclust.rds", "", names(subclust_list))



annotate_clusters <- function(umap_object, subcluster_name, new_cluster_names) {
    # アイデンティティを更新する
    umap_object <- RenameIdents(umap_object, new_cluster_names)

    # meta dataにcell_typeとしてアイデンティティを追加
    umap_object <- AddMetaData(
        umap_object,
        metadata = Idents(umap_object),
        col.name = "celltype"
    )

    ### UMAPをプロット
    DimPlot(
        umap_object,
        reduction = "umap_harmony",
        group.by = "celltype",
        label = TRUE, pt.size = 0.5
    )
    ggsave("output/plot/umap/umap_celltype.png", width = 10, height = 10, dpi = 300)

    ### umap_objectを保存
    saveRDS(umap_object, paste0("data/processed/subclusters/", subcluster_name, "_subclust.rds"))
}

validate_markers <- function(umap_object, subcluster_name, cluster_specific_markers) {
    ### DotPlotでAnnotationが妥当か確認
    # 各クラスターのDotPlotを作成
    for (cluster in names(cluster_specific_markers)) {
        features <- cluster_specific_markers[[cluster]]

        # umap_objectのSCT assayに存在する遺伝子のみをフィルタリング
        valid_features <- features[features %in% rownames(umap_object[["SCT"]]@data)]

        # 有効な遺伝子がない場合はスキップ
        if (length(valid_features) == 0) {
            warning(paste("Cluster", cluster, "には有効な遺伝子がありません。"))
            next
        }

        plot <- DotPlot(umap_object, features = valid_features) +
            ggtitle(paste("Cluster", cluster)) +
            theme_classic()
        ggsave(
            paste0("output/plot/dotplot/subclusters_subcluster/dotplot_cluster_", subcluster_name, "_", cluster, "_markers.png"), 
            plot = plot, width = 10, height = 10, dpi = 300
        )
    }
}

# 関数を呼び出す例
new_cluster_names <- c(
    "1" = "T_cells",
    "2" = "Epithelial_cells",
    "3" = "B_cells",
    "4" = "Plasma_cells",
    "5" = "Endothelial_cells",
    "6" = "Macrophages",
    "7" = "Mast_cells"
)
names(new_cluster_names) <- levels(umap_object)

cluster_specific_markers <- list(
    "T cells" = c("CD3D", "CD3E", "CD4", "CD8A", "IL2RA", "CD27", "CD28", "CTLA4", "PDCD1", "IL2"),
    "Epithelial cells" = c("EPCAM", "KRT18", "KRT19", "CDH1", "MUC1", "KRT7", "KRT14", "KRT8", "CLDN4", "TP63"),
    "B cells" = c("CD19", "CD20", "MS4A1", "CD22", "CD79A", "CD24", "CD38", "IGM", "IGD", "CD27"),
    "Plasma cells" = c("IGHA1", "IGKC", "CD27", "CD38", "XBP1", "PRDM1", "MZB1", "SDC1", "IGHG1", "CD138"),
    "Endothelial cells" = c("PECAM1", "VWF", "CD34", "CD31", "ENG", "CD144", "KDR", "FLT1", "THBD", "ESM1"),
    "Macrophages" = c("CD68", "CD163", "CSF1R", "IL1B", "TNF", "CD14", "CD80", "CD86", "MRC1", "IL6"),
    "Mast cells" = c("FCER1A", "TPSB2", "KIT", "CPA3", "IL9", "TPSAB1", "KITLG", "CMA1", "MCT", "IL4")
)

# subclust_listの一つ一つに対して手動でアノテーションを行う
umap_object <- subclust_list[[1]]
annotate_clusters(umap_object, new_cluster_names)
validate_markers(umap_object, cluster_specific_markers)
