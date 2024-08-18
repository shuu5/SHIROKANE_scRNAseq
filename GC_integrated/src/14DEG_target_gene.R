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
tar_genes <- c("ARID5A", "ARID5B", "KRT13", "TGM2", "AXL", "IL6", "MKI67")

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

### DEG analysis (groupのTumorの中で, 各celltypeに対して, 特定のtar_geneの平均値で２グループに分けて)
# DEG解析の関数を定義
perform_deg_analysis <- function(seurat_obj, cell_type, test = "MAST", tar_genes = tar_genes) {
  print(paste("Analyzing", cell_type))
  
  # 指定されたcell_typeをlistから選択
  subset_obj <- subclust_list[[cell_type]]
  
  # Tumorの細胞のみをフィルタリング
  tumor_cells <- subset_obj[, subset_obj$group == "Tumor"]
  
  # 空のデータフレームチェック
  if (ncol(tumor_cells) == 0) {
    print(paste("Warning: Tumor group has no cells for", cell_type))
    return(NULL)
  }
  
  print(paste("Number of cells in Tumor group:", ncol(tumor_cells)))
  
  # Normal vs Tumorの比較を行う前にPrepSCTFindMarkersを実行
  # SCTアッセイの準備
  tumor_cells <- PrepSCTFindMarkers(tumor_cells) 

  # 各tar_geneに対してDEG解析を実行
  deg_results_list <- lapply(tar_genes, function(tar_gene) {
    # tar_geneがデータに存在するか確認
    if (!tar_gene %in% rownames(tumor_cells@assays$SCT$data)) {
      print(paste("tar_gene:", tar_gene, "はデータに存在しません。スキップします。"))
      return(NULL)  # 存在しない場合はNULLを返す
    }
    
    # tar_geneの中央値でhigh群とlow群に分ける
    mean_value <- mean(tumor_cells@assays$SCT$data[tar_gene,])
    
    # mean_valueが0の場合はスキップ
    if (mean_value == 0) {
      print(paste("tar_gene:", tar_gene, "のmean_valueが0のため、グループを分けられません。スキップします。"))
      return(NULL)
    }
    
    tumor_cells$group <- ifelse(tumor_cells@assays$SCT$data[tar_gene,] > mean_value, "High", "Low")
    print(paste0("tar_gene: ", tar_gene))
    print(table(tumor_cells$group))
    print(paste0("mean_value: ", mean_value))

    # 各グループの細胞数を確認
    if (sum(tumor_cells$group == "High") < 10 || sum(tumor_cells$group == "Low") < 10) {
      print(paste("tar_gene:", tar_gene, "のいずれかのグループが10細胞未満のため、スキップします。"))
      return(NULL)
    }

    # High vs Lowの比較を行う
    deg_results <- FindMarkers(tumor_cells, 
                               ident.1 = "High", 
                               ident.2 = "Low", 
                               group.by = "group", 
                               test.use = test,
                               min.pct = 0,
                               logfc.threshold = 0.25,
                               latent.vars = "nCount_SCT"
                               )
    
    # 結果を保存
    tryCatch({
      write.csv(deg_results, file = paste0("output/csv/DEG/subclusters/tar_genes/", cell_type, "_", tar_gene, "_High_vs_Low_", test, ".csv"))
    }, error = function(e) {
      print(paste("Error saving file for", cell_type, ":", e$message))
    })
    
    return(deg_results)
  })

  # deg_results_listにtar_genesの名前をつける
  names(deg_results_list) <- tar_genes

  # 使用したtest名をdeg_results_listに追加
  deg_results_list$test <- test

  return(deg_results_list)
}


# subclust_list内のすべてのcelltypeに対してDEG解析を実行
select_test <- "MAST"
all_deg_results <- lapply(names(subclust_list), function(cell_type) perform_deg_analysis(
    subclust_list[[cell_type]], cell_type, 
    test = select_test, tar_genes = tar_genes
    ))
names(all_deg_results) <- names(subclust_list)





### Volcano plot
# Volcano plotを作成する関数（p_val_adjを使用し、重要な遺伝子にラベルを付ける）
create_volcano_plot <- function(deg_results, cell_type, tar_gene, top_n_genes = 20, fc_threshold = 1, padj_threshold = 0.05, select_label = NULL, test.use = "wilcox") {  
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
    slice_max(order_by = avg_log2FC, n = top_n_genes)
  
  down_genes <- deg_results %>%
    filter(category == "Down") %>%
    slice_max(order_by = -avg_log2FC, n = top_n_genes)
  
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
  write.csv(up_genes, file = paste0("output/csv/DEG/subclusters/tar_genes/", cell_type, "_", tar_gene, "_Up_", test.use, ".csv"))
  write.csv(down_genes, file = paste0("output/csv/DEG/subclusters/tar_genes/", cell_type, "_", tar_gene, "_Down_", test.use, ".csv"))
  
  # プロットを保存
  ggsave(paste0("output/plot/volcano/subclusters/tar_genes/volcano_plot_", cell_type, "_", tar_gene, "_", test.use, ".png"), p, width = 10, height = 10)
  
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
    ggsave(paste0("output/plot/volcano/subclusters/tar_genes/volcano_plot_selected_", cell_type, "_", tar_gene, "_", test.use, ".png"), p_select, width = 10, height = 10)

    # select_labelのp_valなどの値をcsvに保存し返す。
    write.csv(label_genes_select, file = paste0("output/csv/DEG/subclusters/tar_genes/", cell_type, "_", tar_gene, "_", test.use, "_selected.csv"))
  }
}

# すべてのcelltypeのVolcano plotを作成(top_n_genes = 10, fc_threshold = 1, padj_threshold = 0.05, select_label = tar_genes)
lapply(names(all_deg_results), function(cell_type) {
  lapply(tar_genes, function(tar_gene) {
    print(paste0("cell_type:", cell_type, ", tar_gene:", tar_gene, " volcano plot processing..."))
    
    # NULLチェックを追加
    if (is.null(all_deg_results[[cell_type]][[tar_gene]])) {
      print(paste("Warning: ", cell_type, "'s ", tar_gene, " is NULL. Skipping this processing."))
      return(NULL)  # NULLの場合は処理をスキップ
    }
    
    # label_genes_select_listに格納してall_deg_resultsと同じ名前を各要素につける
    create_volcano_plot(
          deg_results = all_deg_results[[cell_type]][[tar_gene]], 
          cell_type = cell_type, 
          tar_gene = tar_gene,
          top_n_genes = 10, 
          fc_threshold = 1, 
          padj_threshold = 0.05,
          select_label = tar_genes,
          test.use = select_test
      )
  })
})






### 各celltypeの上位遺伝子をDot, Feature, VlnPlotでプロットする
# それそれのサブクラスターの上位の遺伝子をcsvから読み込み
plot_genes_list <- list()
for (cell_type in names(all_deg_results)) {
    for (tar_gene in tar_genes) {
        up_genes_file <- paste0("output/csv/DEG/subclusters/tar_genes/", cell_type, "_", tar_gene, "_Up_", test.use, ".csv")
        down_genes_file <- paste0("output/csv/DEG/subclusters/tar_genes/", cell_type, "_", tar_gene, "_Down_", test.use, ".csv")
        
        # ファイルが存在するか確認
        if (!file.exists(up_genes_file) || !file.exists(down_genes_file)) {
            print(paste("Warning: ", cell_type, "の", tar_gene, "に対応するファイルが見つかりません。処理をスキップします。"))
            next  # ファイルが存在しない場合は次のtar_geneへ
        }
        
        up_genes <- read.csv(up_genes_file)
        down_genes <- read.csv(down_genes_file)
        names(up_genes)[1] <- "gene"
        names(down_genes)[1] <- "gene"

        # 結合
        up_down_genes <- rbind(up_genes, down_genes)
        plot_genes_list[[cell_type]][[tar_gene]] <- up_down_genes
    }
}






# 各celltypeごとにプロットを作成
for (cell_type in names(all_deg_results)) {
    for (tar_gene in tar_genes) {
        plot_genes <- plot_genes_list[[cell_type]][[tar_gene]]
    
        ### Dotplot
        print(paste("Dotplotを作成します:", cell_type, "_", tar_gene))
    
        # groupで分ける
        DotPlot(
            subclust_list[[cell_type]],
            features = plot_genes$gene,
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
        ggsave(paste0("output/plot/dotplot/subclusters/tar_genes/dotplot_", cell_type, "_", tar_gene, "_group_", select_test, ".png"), width = 10, height = 10)

        # tissueで分ける
        DotPlot(
            subclust_list[[cell_type]],
            features = plot_genes$gene,
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
        ggsave(paste0("output/plot/dotplot/subclusters/tar_genes/dotplot_", cell_type, "_", tar_gene, "_tissue_", select_test, ".png"), width = 10, height = 10)

        ### FeaturePlot, VlnPlot
        for (gene in plot_genes) {
            # 遺伝子がデータに存在するか確認
            if (gene %in% rownames(subclust_list[[cell_type]]@assays$SCT$data)) {
                print(paste("遺伝子が見つかりました:", gene))

                ### FeaturePlot
                print(paste("FeaturePlotを作成します:", gene))
                FeaturePlot(
                    subclust_list[[cell_type]],
                    features = gene
                )
                ggsave(paste0("output/plot/feature/subclusters/tar_genes/featureplot_", cell_type, "_", tar_gene, "_", gene, "_", select_test, ".png"), width = 15, height = 10)

                # groupでsplit
                FeaturePlot(
                    subclust_list[[cell_type]],
                    features = gene,
                    split.by = "group"
                )
                ggsave(paste0("output/plot/feature/subclusters/tar_genes/featureplot_", cell_type, "_", tar_gene, "_", gene, "_group_", select_test, ".png"), width = 15, height = 10)

                # tissueでsplit
                FeaturePlot(
                    subclust_list[[cell_type]],
                    features = gene,
                    split.by = "tissue"
                )
                ggsave(paste0("output/plot/feature/subclusters/tar_genes/featureplot_", cell_type, "_", tar_gene, "_", gene, "_tissue_", select_test, ".png"), width = 15, height = 10)

                ### Violin plot
                print(paste("Violin plotを作成します:", gene))
                VlnPlot(
                    subclust_list[[cell_type]],
                    features = gene
                )
                ggsave(paste0("output/plot/violin/subclusters/tar_genes/violinplot_", cell_type, "_", tar_gene, "_", gene, "_", select_test, ".png"), width = 15, height = 10)

                # groupでsplit
                VlnPlot(
                    subclust_list[[cell_type]],
                    features = gene,
                    split.by = "group"
                )
                ggsave(paste0("output/plot/violin/subclusters/tar_genes/violinplot_", cell_type, "_", tar_gene, "_", gene, "_group_", select_test, ".png"), width = 15, height = 10)

                # tissueでsplit
                VlnPlot(
                    subclust_list[[cell_type]],
                    features = gene,
                    split.by = "tissue"
                )
                ggsave(paste0("output/plot/violin/subclusters/tar_genes/violinplot_", cell_type, "_", tar_gene, "_", gene, "_tissue_", select_test, ".png"), width = 15, height = 10)
            } else {
            warning(paste("遺伝子が見つかりません:", gene))
        }
    }
    }
}