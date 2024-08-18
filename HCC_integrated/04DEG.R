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

# 各meta.dataの情報を確認
merged_object
head(merged_object)
table(merged_object@meta.data$tissue)
table(merged_object@meta.data$split_param)
table(merged_object@meta.data$etiology)
table(merged_object@meta.data$cell_type)
table(merged_object@meta.data$ordered_cell_type)

# 結果保存用のディレクトリを作成
dir.create(paste0("results/", folder_name, "/DEG"), recursive = TRUE, showWarnings = FALSE)


subset_obj <- subset(merged_object, cell_type == "T_cell_1")

subset_obj <- PrepSCTFindMarkers(subset_obj)

deg_results <- FindMarkers(subset_obj, 
                             ident.1 = "Tumor", 
                             ident.2 = "Normal", 
                             group.by = "tissue", 
                             test.use = "wilcox",
                             min.pct = 0.1,
                             logfc.threshold = 0.25)

# DEG解析の関数を定義
perform_deg_analysis <- function(seurat_obj, ct) {
  print(paste("Analyzing", ct))
  
  # 指定されたcell_typeのみを抽出
  subset_obj <- subset(seurat_obj, cell_type == ct)
  
  print(paste("Number of cells in", ct, ":", ncol(subset_obj)))
  print(table(subset_obj$tissue))
  
  # Normal vs Tumorの比較を行う前にPrepSCTFindMarkersを実行
  # SCTアッセイの準備
  subset_obj <- PrepSCTFindMarkers(subset_obj)  # latent.varsを削除

  # Normal vs Tumorの比較を行う
  deg_results <- FindMarkers(subset_obj, 
                             ident.1 = "Tumor", 
                             ident.2 = "Normal", 
                             group.by = "tissue", 
                             test.use = "wilcox",
                             min.pct = 0.1,
                             logfc.threshold = 0.25,
                             latent.vars = "nCount_SCT")
  
  # 結果を保存
  tryCatch({
    write.csv(deg_results, file = paste0("results/", folder_name, "/DEG/DEG_", ct, "_Tumor_vs_Normal.csv"))
  }, error = function(e) {
    print(paste("Error saving file for", ct, ":", e$message))
  })
  
  return(deg_results)
}

# cell_typesを取得
cell_types <- unique(merged_object@meta.data$cell_type)
print("Cell types to analyze:")
print(cell_types)

# 各cell_typeに対してDEG解析を実行
all_deg_results <- list()
for (ct in cell_types) {
  print(paste("Starting analysis for", ct))
  all_deg_results[[ct]] <- perform_deg_analysis(merged_object, ct)
  print(paste("Finished analysis for", ct))
}

# 結果のリストに名前を付ける
names(all_deg_results) <- cell_types

# 結果の確認（例：上位10個の遺伝子を表示）
lapply(names(all_deg_results), function(ct) {
  print(paste("Top 10 DEGs for", ct))
  print(head(all_deg_results[[ct]], 10))
})









### plot
# Volcano plotを作成する関数（p_val_adjを使用し、重要な遺伝子にラベルを付ける）
create_volcano_plot <- function(deg_results, ct, top_n_genes = 20, fc_threshold = 1, padj_threshold = 0.05) {
  # NAを除去し、adj_pが0でない遺伝子のみを保持
  deg_results <- deg_results %>%
    na.omit() %>%
    filter(p_val_adj > 0)
  
  # log10(p_val_adj)を計算
  deg_results$log10_padj <- -log10(deg_results$p_val_adj)
  
  # 各遺伝子にカテゴリを割り当て
  deg_results$category <- case_when(
    deg_results$p_val_adj < padj_threshold & deg_results$avg_log2FC > fc_threshold ~ "Up",
    deg_results$p_val_adj < padj_threshold & deg_results$avg_log2FC < -fc_threshold ~ "Down",
    TRUE ~ "NS"
  )
  
  # 発現が上昇している遺伝子と減少している遺伝子を別々に選択
  up_genes <- deg_results %>%
    filter(category == "Up") %>%
    top_n(top_n_genes, wt = avg_log2FC)
  
  down_genes <- deg_results %>%
    filter(category == "Down") %>%
    top_n(top_n_genes, wt = -avg_log2FC)
  
  # 選択した遺伝子を結合
  label_genes <- rbind(up_genes, down_genes)
  
  # プロットの作成
  p <- ggplot(deg_results, aes(x = avg_log2FC, y = log10_padj, color = category)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    theme_classic() +
    theme(legend.position = "right") +
    labs(title = paste("Volcano Plot -", ct),
         x = "log2 Fold Change",
         y = "-log10 adjusted p-value") +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed") +
    geom_text_repel(
      data = label_genes,
      aes(label = rownames(label_genes)),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      point.padding = 0.5,
      segment.color = "grey50"
    )
  
  # プロットを保存
  ggsave(paste0("results/", folder_name, "/DEG/volcano_plot_", ct, ".png"), p, width = 12, height = 10)
  
  return(list(plot = p, labeled_genes = rownames(label_genes)))
}

# Feature plotを作成する関数
create_feature_plots <- function(seurat_obj, genes, ct) {
  print(paste("Creating Feature plots for", ct))
  print(paste("Number of genes to plot:", length(genes)))
  
  # 指定されたcell_typeのみを抽出
  subset_obj <- subset(seurat_obj, ordered_cell_type == ct)
  
  # Feature plotを作成
  p <- FeaturePlot(subset_obj, 
                   features = genes, 
                   split.by = "tissue",
                   ncol = 4,  # 1行あたりの遺伝子数を増やす
                   reduction = "umap")
  
  # プロットを保存
  ggsave(paste0("results/", folder_name, "/DEG/feature_plot_", ct, ".png"), 
         p, 
         width = 20,  # 幅を調整
         height = 10 * ceiling(length(genes)),  # 高さを調整
         limitsize = FALSE)  # サイズ制限を解除
  
  print(paste("Feature plots for", ct, "have been saved."))
  
  return(p)
}

# Volcano plotとFeature plotを作成
volcano_and_feature_plots <- lapply(names(all_deg_results), function(ct) {
  print(paste("Processing", ct))
  
  # Volcano plotを作成
  volcano_result <- create_volcano_plot(all_deg_results[[ct]], ct, top_n_genes = 5)
  
  # ラベル付けされた遺伝子を取得
  labeled_genes <- volcano_result$labeled_genes
  
  print(paste("Number of labeled genes for", ct, ":", length(labeled_genes)))
  
  # Feature plotを作成
  feature_plot <- create_feature_plots(merged_object, labeled_genes, ct)
  
  return(list(volcano_plot = volcano_result$plot, feature_plot = feature_plot))
})

# 結果のリストに名前を付ける
names(volcano_and_feature_plots) <- cell_types

print("Volcano plots and Feature plots have been created for all cell types.")