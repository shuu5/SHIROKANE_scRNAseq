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
library(MAST)
library(ggrepel)

# folder name
folder_name <- "HCC_ALL_integrated"

# harmony_merged_objectを読み込み
merged_object <- readRDS(paste0("data/", folder_name, "/merged_object.rds"))

# tissueごとにHepatic_stellate_cellのデータをフィルタリング
hepatocyte_data <- subset(merged_object, cell_type == "Hepatic_Stellate_Cell")

# NormalとTumorでの細胞タイプを分ける
hepatocyte_data$celltype.tissue <- paste(hepatocyte_data$cell_type, hepatocyte_data$tissue, sep = "_")

# データを準備
hepatocyte_data <- PrepSCTFindMarkers(hepatocyte_data)  # データを準備

# NormalとTumorの間での差異を調べる
Idents(hepatocyte_data) <- "celltype.tissue"
differential_genes <- FindMarkers(hepatocyte_data, ident.1 = "Hepatic_Stellate_Cell_Normal", ident.2 = "Hepatic_Stellate_Cell_Tumor", verbose = FALSE)

# 結果を表示
head(differential_genes, n = 10)

# NormalとTumorの間での遺伝子発現を集約
aggregated_data <- AggregateExpression(hepatocyte_data, 
                                       group.by = c("cell_type", "tissue"), 
                                       return.seurat = TRUE)

# TumorとNormalの発現差を計算
expression_diff <- aggregated_data[["SCT"]]@data[, "Hepatic-Stellate-Cell_Tumor"] - aggregated_data[["SCT"]]@data[, "Hepatic-Stellate-Cell_Normal"]

# Tumor - Normalの値が最も大きい上位10遺伝子を選択
top_positive_genes <- head(rownames(aggregated_data)[order(expression_diff, decreasing = TRUE)], 10)

# Tumor - Normalの値が最も小さい上位10遺伝子を選択
top_negative_genes <- head(rownames(aggregated_data)[order(expression_diff)], 10)

# Hepatic Stellate Cellの散布図を作成
p1 <- ggplot(data = as.data.frame(aggregated_data[["SCT"]]@data), aes(x = `Hepatic-Stellate-Cell_Normal`, y = `Hepatic-Stellate-Cell_Tumor`)) +
      geom_point(aes(color = ifelse(rownames(aggregated_data) %in% top_positive_genes, "red", 
                                     ifelse(rownames(aggregated_data) %in% top_negative_genes, "blue", "grey"))), 
                 size = 1) +
      scale_color_identity() +  # 色をそのまま使用
      theme_classic()  # classicテーマを適用

# 赤い点（差が大きい方）と青い点（差が小さい方）をラベル付け
p2 <- p1 + 
      ggrepel::geom_text_repel(data = subset(as.data.frame(aggregated_data[["SCT"]]@data), 
                                              rownames(aggregated_data) %in% top_positive_genes), 
                                aes(label = rownames(aggregated_data)[rownames(aggregated_data) %in% top_positive_genes]), 
                                color = "red", 
                                size = 3, 
                                nudge_x = 0, 
                                nudge_y = 0)

p3 <- p2 + 
      ggrepel::geom_text_repel(data = subset(as.data.frame(aggregated_data[["SCT"]]@data), 
                                              rownames(aggregated_data) %in% top_negative_genes), 
                                aes(label = rownames(aggregated_data)[rownames(aggregated_data) %in% top_negative_genes]), 
                                color = "blue", 
                                size = 3, 
                                nudge_x = 0, 
                                nudge_y = 0)

# プロットを表示
print(p3)

# プロットを保存
ggsave(
    paste0("plot/", folder_name, "/aggregate_tumor_normal/Hepatic_Stellate_Cell_scatter.png"), 
    plot = p3, width = 10, height = 10, dpi = 300
    )

# top_positive_genesとtop_negative_genesをcsvに保存
write.csv(top_positive_genes, paste0("data/", folder_name, "/top_positive_genes.csv"), row.names = FALSE)
write.csv(top_negative_genes, paste0("data/", folder_name, "/top_negative_genes.csv"), row.names = FALSE)



# Positive遺伝子のフィーチャープロット（tissueで分割）
p_positive <- FeaturePlot(merged_object, features = top_positive_genes, 
                           cols = c("grey", "red"), 
                           ncol = 2, 
                           split.by = "tissue") + 
              theme_classic() + 
              ggtitle("Top Positive Genes")

# Negative遺伝子のフィーチャープロット（tissueで分割）
p_negative <- FeaturePlot(merged_object, features = top_negative_genes, 
                           cols = c("grey", "blue"), 
                           ncol = 2, 
                           split.by = "tissue") + 
              theme_classic() + 
              ggtitle("Top Negative Genes")

# プロットを表示
print(p_positive)
print(p_negative)

# プロットを保存
ggsave(
    paste0("plot/", folder_name, "/aggregate_tumor_normal/Hepatic_Stellate_Cell_feature_plot_positive.png"), 
    plot = p_positive, width = 10, height = 30, dpi = 300
    )

ggsave(
    paste0("plot/", folder_name, "/aggregate_tumor_normal/Hepatic_Stellate_Cell_feature_plot_negative.png"), 
    plot = p_negative, width = 10, height = 30, dpi = 300
    )



# Positive遺伝子のバイオリンプロット（tissueで分割）
v_positive <- VlnPlot(merged_object, features = top_positive_genes, 
                      group.by = "cell_type",
                      split.by = "tissue",
                      cols = c("grey", "red"), 
                      ncol = 2) + 
              theme_classic() + 
              ggtitle("Top Positive Genes Vln Plot")

# Negative遺伝子のバイオリンプロット（tissueで分割）
v_negative <- VlnPlot(merged_object, features = top_negative_genes, 
                      group.by = "cell_type",
                      split.by = "tissue",
                      cols = c("grey", "blue"), 
                      ncol = 2) + 
              theme_classic() + 
              ggtitle("Top Negative Genes Vln Plot")

# プロットを表示
print(v_positive)
print(v_negative)

# プロットを保存
ggsave(
    paste0("plot/", folder_name, "/aggregate_tumor_normal/Hepatic_Stellate_Cell_vln_plot_positive.png"), 
    plot = v_positive, width = 10, height = 30, dpi = 300
    )

ggsave(
    paste0("plot/", folder_name, "/aggregate_tumor_normal/Hepatic_Stellate_Cell_vln_plot_negative.png"), 
    plot = v_negative, width = 10, height = 30, dpi = 300
    )
