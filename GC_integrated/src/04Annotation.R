library(Seurat)
library(tidyverse)

# umap_objectを読み込み
umap_object <- readRDS("data/processed/umap_object.rds")

head(umap_object@meta.data)

# Identを設定
Idents(umap_object) <- umap_object@meta.data$SCT_snn_res.0.1

# レベルを並び替え
umap_object <- SetIdent(umap_object, value = factor(Idents(umap_object), levels = as.character(1:15)))

### 各クラスタの変動遺伝子からクラスターに名前を付ける
new_cluster_names <- c(
    "1" = "T_cells",
    "2" = "B_cells",
    "3" = "T_cells",
    "4" = "Epithelial_cells",
    "5" = "Plasma_cells",
    "6" = "Epithelial_cells",
    "7" = "Plasma_cells",
    "8" = "Macrophages",
    "9" = "B_cells",
    "10" = "Smooth_muscle_cells",
    "11" = "Mast_cells",
    "12" = "Endothelial_cells",
    "13" = "Plasma_cells",
    "14" = "Epithelial_cells",
    "15" = "Plasma_cells"
)
names(new_cluster_names) <- levels(umap_object)

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



### DotPlotでAnnotationが妥当か確認
# 各クラスターの特異的なマーカー遺伝子
cluster_specific_markers <- list(
    "T cells" = c("CD3D", "CD3E", "CD4", "CD8A", "IL2RA", "CD27", "CD28", "CTLA4", "PDCD1", "IL2"),
    "Epithelial cells" = c("EPCAM", "KRT18", "KRT19", "CDH1", "MUC1", "KRT7", "KRT14", "KRT8", "CLDN4", "TP63"),
    "B cells" = c("CD19", "CD20", "MS4A1", "CD22", "CD79A", "CD24", "CD38", "IGM", "IGD", "CD27"),
    "Plasma cells" = c("IGHA1", "IGKC", "CD27", "CD38", "XBP1", "PRDM1", "MZB1", "SDC1", "IGHG1", "CD138"),
    "Endothelial cells" = c("PECAM1", "VWF", "CD34", "CD31", "ENG", "CD144", "KDR", "FLT1", "THBD", "ESM1"),
    "Macrophages" = c("CD68", "CD163", "CSF1R", "IL1B", "TNF", "CD14", "CD80", "CD86", "MRC1", "IL6"),
    "Mast cells" = c("FCER1A", "TPSB2", "KIT", "CPA3", "IL9", "TPSAB1", "KITLG", "CMA1", "MCT", "IL4"), 
    "Smooth muscle cells" = c("ACTA2", "MYH11", "CNN1", "TAGLN", "SM22", "MYLK", "CALD1", "DES", "VIM", "FAP"), 
    "NK cells" = c("NKG7", "CD56", "KLRD1", "GNLY", "PRF1", "GZMB", "CD16", "IL2RB", "NKG2A", "FCGR3A"),
    "Monocytes" = c("CD14", "CD16", "CD68", "CSF1R", "CCR2", "IL1B", "TNF", "CD163", "CD86", "MRC1")
)

# 各クラスターのDotPlotを作成
for (cluster in names(cluster_specific_markers)) {
    print(paste0("Cluster: ", cluster))

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
        paste0("output/plot/dotplot/dotplot_cluster_", cluster, "_markers.png"), 
        plot = plot, width = 10, height = 10, dpi = 300
    )
}

# 各クラスターのFeaturePlotを作成
for (cluster in names(cluster_specific_markers)) {
    print(paste0("Cluster: ", cluster))
    features <- cluster_specific_markers[[cluster]]
    plot <- FeaturePlot(umap_object, features = features) +
        ggtitle(paste("Cluster", cluster)) +
        theme_classic()
    ggsave(
        paste0("output/plot/feature/validation/feature_cluster_", cluster, "_markers.png"), 
        plot = plot, width = 15, height = 10, dpi = 300
    )
}


### umap_objectを保存
saveRDS(umap_object, "data/processed/umap_object.rds")

