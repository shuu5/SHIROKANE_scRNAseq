library(Seurat)
library(tidyverse)
library(presto)
library(ggrepel)
library(MAST)


# data/processed/subclusters/内のsubclusteringしたオブジェクトを読み込み
subclust_list <- lapply(list.files("data/processed/subclusters/", pattern = "*.rds", full.names = TRUE), readRDS)
names(subclust_list) <- gsub("data/processed/subclusters/(.*)_subclust.rds", "\\1", list.files("data/processed/subclusters/", pattern = "*.rds"))
names(subclust_list) <- gsub("_subclust.rds", "", names(subclust_list))

# 特に見たいgenes
tar_genes <- c("ARID5A", "ARID5B", "KRT13", "TGM2", "AXL", "MKI67")


# 各subclust_listのオブジェクトにtar_genesが含まれているか確認
lapply(names(subclust_list), function(cell_type) {
  subset_obj <- subclust_list[[cell_type]]
  missing_genes <- setdiff(tar_genes, rownames(subset_obj))
  
  if (length(missing_genes) > 0) {
    print(paste(cell_type, "には以下の遺伝子が含まれていません:", paste(missing_genes, collapse = ", ")))
  } else {
    print(paste(cell_type, "にはすべてのtar_genesが含まれています。"))
  }
})


### DEG analysis (Tumor vs Paratumor, 各celltypeに対して)
# DEG解析の関数を定義
perform_deg_analysis <- function(seurat_obj, cell_type, test = "wilcox") {
  print(paste("Analyzing", cell_type))
  
  # 指定されたcell_typeをlistから選択
  subset_obj <- subclust_list[[cell_type]]
  
  print(paste("Number of cells in", cell_type, ":", ncol(subset_obj)))
  print(table(subset_obj$group))
  
  # Normal vs Tumorの比較を行う前にPrepSCTFindMarkersを実行
  # SCTアッセイの準備
  subset_obj <- PrepSCTFindMarkers(subset_obj) 

  # Normal vs Tumorの比較を行う
  deg_results <- FindMarkers(subset_obj, 
                             ident.1 = "Tumor", 
                             ident.2 = "Non_Tumor", 
                             group.by = "group", 
                             test.use = test,
                             min.pct = 0,
                             logfc.threshold = 0.25,
                             latent.vars = "nCount_SCT"
                             )
  
  # 結果を保存
  tryCatch({
    write.csv(deg_results, file = paste0("output/csv/DEG/subclusters/", cell_type, "_Tumor_vs_Non-Tumor_", test, ".csv"))
  }, error = function(e) {
    print(paste("Error saving file for", cell_type, ":", e$message))
  })
  
  return(deg_results)
}

# subclust_list内のすべてのcelltypeに対してDEG解析を実行
select_test <- "MAST"
all_deg_results <- lapply(names(subclust_list), function(cell_type) perform_deg_analysis(subclust_list[[cell_type]], cell_type, test = select_test))
names(all_deg_results) <- names(subclust_list)

# 使用したtest名をall_deg_resultsに追加
all_deg_results <- lapply(all_deg_results, function(deg_result) {
  deg_result$test <- select_test
  return(deg_result)
})



### Volcano plot
# Volcano plotを作成する関数（p_val_adjを使用し、重要な遺伝子にラベルを付ける）
create_volcano_plot <- function(deg_results, cell_type, top_n_genes = 20, fc_threshold = 1, padj_threshold = 0.05, select_label = NULL, test.use = "wilcox") {  
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
    labs(title = paste("Volcano Plot -", cell_type),
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

  # 発現が上昇・減少している遺伝子をcsvに保存
  write.csv(up_genes, file = paste0("output/csv/DEG/subclusters/", cell_type, "_Up_", test.use, ".csv"))
  write.csv(down_genes, file = paste0("output/csv/DEG/subclusters/", cell_type, "_Down_", test.use, ".csv"))
  
  # プロットを保存
  ggsave(paste0("output/plot/volcano/subclusters/volcano_plot_", cell_type, "_", test.use, ".png"), p, width = 12, height = 10)
  
  # select_labelが指定されている場合、追加のプロットを作成
  if (!is.null(select_label)) {
    # select_labelに指定された遺伝子の存在を確認
    missing_genes <- setdiff(select_label, rownames(deg_results))
    if (length(missing_genes) > 0) {
      print(paste("以下の遺伝子はデータに存在しません:", paste(missing_genes, collapse = ", ")))
    }
    
    label_genes_select <- deg_results[rownames(deg_results) %in% select_label, ]
    
    p_select <- ggplot(deg_results, aes(x = avg_log2FC, y = log10_padj, color = category)) +
      geom_point(alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Up" = "red", "Down" = "blue", "label" = "darkgreen", "NS" = "grey")) +
      theme_classic() +
      theme(legend.position = "right") +
      labs(title = paste("Volcano Plot with Selected Genes -", cell_type),
           x = "log2 Fold Change",
           y = "-log10 adjusted p-value") +
      geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed") +
      geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed") +
      geom_point(
        data = label_genes_select,
        aes(x = avg_log2FC, y = log10_padj, color = "label"),
        size = 3,
        alpha = 0.6
      ) +
      geom_text_repel(
        data = label_genes_select,
        aes(label = rownames(label_genes_select), color = "label"),
        size = 3,
        max.overlaps = 20,
        box.padding = 0.5,
        point.padding = 0.5,
        segment.color = "grey50"
      )
    
    # 追加のプロットを保存
    ggsave(paste0("output/plot/volcano/subclusters/volcano_plot_selected_", cell_type, "_", test.use, ".png"), p_select, width = 12, height = 10)

    # select_labelのp_valなどの値をcsvに保存し返す。
    write.csv(label_genes_select, file = paste0("output/csv/DEG/subclusters/", cell_type, "_", test.use, "_selected.csv"))
    return(label_genes_select)
  }
}

# すべてのcelltypeのVolcano plotを作成(top_n_genes = 10, fc_threshold = 1, padj_threshold = 0.05, select_label = tar_genes)
lapply(names(all_deg_results), function(cell_type) {
  # label_genes_select_listに格納してall_deg_resultsと同じ名前を各要素につける
  label_genes_select_list[[cell_type]] <- create_volcano_plot(
        deg_results = all_deg_results[[cell_type]], 
        cell_type = cell_type, 
        top_n_genes = 10, 
        fc_threshold = 1, 
        padj_threshold = 0.05,
        select_label = tar_genes,
        test.use = select_test
    )
})

print(label_genes_select_list)

### 各celltypeの上位遺伝子をDot, Feature, VlnPlotでプロットする
# それそれのサブクラスターの上位の遺伝子をcsvから読み込み
tar_genes_list <- list()
for (cell_type in names(all_deg_results)) {
    up_genes <- read.csv(paste0("output/csv/DEG/subclusters/", cell_type, "_Up_", select_test, ".csv"))
    down_genes <- read.csv(paste0("output/csv/DEG/subclusters/", cell_type, "_Down_", select_test, ".csv"))
    names(up_genes)[1] <- "gene"
    names(down_genes)[1] <- "gene"

    # 結合
    tar_genes <- rbind(up_genes, down_genes)
    tar_genes_list[[cell_type]] <- tar_genes
}



# 各celltypeごとにプロットを作成
for (cell_type in names(all_deg_results)) {
    tar_genes <- tar_genes_list[[cell_type]]$gene
    
    ### Dotplot
    print(paste("Dotplotを作成します:", cell_type))
    # groupで分ける
    DotPlot(
        subclust_list[[cell_type]],
        features = tar_genes,
        group.by = "group"
    ) + 
    coord_flip() +
    theme_classic() +
    theme(
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 12),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
    )
    ggsave(paste0("output/plot/dotplot/subclusters/dotplot_", cell_type, "_group_", select_test, ".png"), width = 10, height = 10)

    # tissueで分ける
    DotPlot(
        subclust_list[[cell_type]],
        features = tar_genes,
        group.by = "tissue"
    ) +
    coord_flip() +
    theme_classic() +
    theme(
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 12),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
    )
    ggsave(paste0("output/plot/dotplot/subclusters/dotplot_", cell_type, "_tissue_", select_test, ".png"), width = 10, height = 10)

    ### FeaturePlot, VlnPlot
    for (gene in tar_genes) {
        # 遺伝子がデータに存在するか確認
        if (gene %in% rownames(subclust_list[[cell_type]]@assays$SCT$data)) {
            print(paste("遺伝子が見つかりました:", gene))

            ### FeaturePlot
            print(paste("FeaturePlotを作成します:", gene))
            FeaturePlot(
                subclust_list[[cell_type]],
                features = gene
            )
            ggsave(paste0("output/plot/feature/subclusters/featureplot_", cell_type, "_", gene, "_", select_test, ".png"), width = 15, height = 10)

            # groupでsplit
            FeaturePlot(
                subclust_list[[cell_type]],
                features = gene,
                split.by = "group"
            )
            ggsave(paste0("output/plot/feature/subclusters/featureplot_", cell_type, "_", gene, "_group_", select_test, ".png"), width = 15, height = 10)

            # tissueでsplit
            FeaturePlot(
                subclust_list[[cell_type]],
                features = gene,
                split.by = "tissue"
            )
            ggsave(paste0("output/plot/feature/subclusters/featureplot_", cell_type, "_", gene, "_tissue_", select_test, ".png"), width = 15, height = 10)

            ### Violin plot
            print(paste("Violin plotを作成します:", gene))
            VlnPlot(
                subclust_list[[cell_type]],
                features = gene
            )
            ggsave(paste0("output/plot/violin/subclusters/violinplot_", cell_type, "_", gene, "_", select_test, ".png"), width = 15, height = 10)

            # groupでsplit
            VlnPlot(
                subclust_list[[cell_type]],
                features = gene,
                split.by = "group"
            )
            ggsave(paste0("output/plot/violin/subclusters/violinplot_", cell_type, "_", gene, "_group_", select_test, ".png"), width = 15, height = 10)

            # tissueでsplit
            VlnPlot(
                subclust_list[[cell_type]],
                features = gene,
                split.by = "tissue"
            )
            ggsave(paste0("output/plot/violin/subclusters/violinplot_", cell_type, "_", gene, "_tissue_", select_test, ".png"), width = 15, height = 10)
        } else {
            warning(paste("遺伝子が見つかりません:", gene))
        }
    }
}