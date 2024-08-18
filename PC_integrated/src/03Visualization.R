library(Seurat)
library(tidyverse)
library(rstatix)

### UMAPの可視化
# celltypeで色分け
DimPlot(seurat_obj_pdac, reduction = "umap_harmony", group.by = "celltype")
ggsave("output/plot/UMAP/UMAP_celltype.png")

# celltypeで色分け、groupでsplit
DimPlot(seurat_obj_pdac, reduction = "umap_harmony", group.by = "celltype", split.by = "group")
ggsave("output/plot/UMAP/UMAP_celltype_group.png", width = 15, height = 10)



### 単遺伝子の可視化
FeaturePlot(seurat_obj_pdac, features = c("IL33"), split.by = "group")
ggsave("output/plot/Feature/FeaturePlot_IL33_group.png", width = 10, height = 10)

VlnPlot(seurat_obj_pdac, features = c("IL33"), split.by = "group", group.by = "celltype")
ggsave("output/plot/Vln/VlnPlot_IL33_group_celltype.png", width = 10, height = 10)

RidgePlot(seurat_obj_pdac, features = "IL33", group.by = "celltype", slot = "data")
ggsave("output/plot/Ridge/RidgePlot_IL33.png", width = 10, height = 8)



### 複数遺伝子の可視化
tar_feature <- c("IL33", "ARID5A", "ARID5B")

## DotPlot
DotPlot(
    seurat_obj_pdac, features = tar_feature, 
    group.by = "celltype", split.by = "group", 
    cols = c("blue", "blue")
    )
ggsave("output/plot/Dot/DotPlot_IL33.png", width = 10, height = 8)




## PseudoBulkによる可視化とDEG解析
library(Seurat)
library(tidyverse)

# 見たい遺伝子（特徴）を指定
tar_feature <- c("IL33", "ARID5A")

# 結果を格納するリストを初期化
test_results <- list()

# 各セルタイプに対して検定を実行
for (cell_type in unique(seurat_obj_pdac$celltype)) {
  # セルタイプでサブセット
  subset_obj <- subset(seurat_obj_pdac, celltype == cell_type)
  
  # サブセットに対してPrepSCTFindMarkersを再度実行
  subset_obj <- PrepSCTFindMarkers(subset_obj)
  
  # グループ間で差異のある遺伝子を検出
  markers <- FindMarkers(subset_obj, 
                         ident.1 = "Tumor", 
                         ident.2 = "Adjacent", 
                         group.by = "group", 
                         features = tar_feature,
                         test.use = "wilcox",
                         slot = "data",
                         assay = "SCT")
  
  # 結果をリストに追加
  test_results[[cell_type]] <- markers
}

# 結果をデータフレームにまとめる
result_df <- do.call(rbind, lapply(names(test_results), function(cell_type) {
  df <- data.frame(celltype = cell_type, feature = rownames(test_results[[cell_type]]), test_results[[cell_type]], row.names = NULL)
  df
}))

# p値で並べ替え
result_df <- result_df[order(result_df$p_val), ]

# 結果を表示
print(result_df)

# 結果をCSVファイルに保存
write.csv(result_df, file = "output/differential_expression_FindMarkers.csv", row.names = FALSE)

# 各遺伝子に対してバープロットを作成
for (feature in tar_feature) {
  feature_results <- result_df %>% filter(feature == !!feature)
  
  p <- ggplot(feature_results, aes(x = celltype, y = -log10(p_val), fill = avg_log2FC)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "細胞タイプ", 
         y = "-log10(p-value)", 
         fill = "平均 log2 発現量変化",
         title = paste(feature, "の細胞タイプ別発現差"))
  
  ggsave(paste0("output/plot/DEG_BarPlot_FindMarkers_", feature, "_by_celltype.png"), plot = p, width = 12, height = 8)
}