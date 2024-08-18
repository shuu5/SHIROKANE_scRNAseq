library(Seurat)
library(tidyverse)
library(presto)
library(ggrepel)

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



### subclust_listの各オブジェクトのgroup間の発現差をpseudo-bulkで計算しplotとcsvで保存
for (i in seq_along(subclust_list)) {
    subclust <- subclust_list[[i]]
    
    # 各サブクラスタの平均発現値を計算
    agg_subclust <- AggregateExpression(subclust, group.by = c("group"), return.seurat = TRUE)

    # 発言の差を計算
    expression_diff <- agg_subclust@assays$SCT$data[, paste0("Tumor")] - agg_subclust@assays$SCT$data[, paste0("Non-Tumor")]

    # 発言の差が大きいものを上位10個
    top_positive_genes <- head(rownames(agg_subclust)[order(expression_diff, decreasing = TRUE)], 10)
    # 発言の差が小さいものを上位10個
    top_negative_genes <- head(rownames(agg_subclust)[order(expression_diff)], 10)


    # 散布図を作成
    p1 <- ggplot(
        data = as.data.frame(agg_subclust@assays$SCT$data), 
        aes(x = agg_subclust@assays$SCT$data[, "Non-Tumor"], y = agg_subclust@assays$SCT$data[, "Tumor"]) 
    ) +
          geom_point(aes(color = ifelse(rownames(agg_subclust) %in% top_positive_genes, "red", 
                                         ifelse(rownames(agg_subclust) %in% top_negative_genes, "blue", "grey"))), 
                     size = 1) +
          labs(x = "Non-Tumor", y = "Tumor") +
          scale_color_identity() +  # 色をそのまま使用
          theme_classic()

    # 赤い点（差が大きい方）と青い点（差が小さい方）をラベル付け
    p2 <- p1 + 
          ggrepel::geom_text_repel(data = subset(as.data.frame(agg_subclust@assays$SCT$data), 
                                                  rownames(agg_subclust) %in% top_positive_genes), 
                                aes(x = agg_subclust@assays$SCT$data[rownames(agg_subclust) %in% top_positive_genes, "Non-Tumor"],
                                    y = agg_subclust@assays$SCT$data[rownames(agg_subclust) %in% top_positive_genes, "Tumor"],
                                    label = rownames(agg_subclust)[rownames(agg_subclust) %in% top_positive_genes]), 
                                color = "red", 
                                size = 3, 
                                nudge_x = 0, 
                                nudge_y = 0)

    p3 <- p2 + 
          ggrepel::geom_text_repel(data = subset(as.data.frame(agg_subclust@assays$SCT$data), 
                                                  rownames(agg_subclust) %in% top_negative_genes), 
                                aes(x = agg_subclust@assays$SCT$data[rownames(agg_subclust) %in% top_negative_genes, "Non-Tumor"],
                                    y = agg_subclust@assays$SCT$data[rownames(agg_subclust) %in% top_negative_genes, "Tumor"],
                                    label = rownames(agg_subclust)[rownames(agg_subclust) %in% top_negative_genes]), 
                                color = "blue", 
                                size = 3, 
                                nudge_x = 0, 
                                nudge_y = 0)

    # プロットを保存
    ggsave(
        paste0("output/plot/scatter/pseudobulk/", names(subclust_list)[i], ".png"), 
        plot = p3, width = 10, height = 10, dpi = 300
    )

    # top_positive_genesとtop_negative_genesをcsvに保存
    write.csv(top_positive_genes, paste0("output/csv/pseudobulk/", names(subclust_list)[i], "_top_positive_genes.csv"), row.names = FALSE)
    write.csv(top_negative_genes, paste0("output/csv/pseudobulk/", names(subclust_list)[i], "_top_negative_genes.csv"), row.names = FALSE)

    # 各Positive遺伝子のフィーチャープロット（groupで分割）
    for (gene in top_positive_genes) {
        print(paste("Positive Gene:", gene, "feature plot"))
        p_positive <- FeaturePlot(subclust, features = gene, 
                                   cols = c("grey", "red"), 
                                   ncol = 1, 
                                   split.by = "group") + 
                      theme_classic() + 
                      ggtitle(paste("Positive Gene:", gene))

        # プロットを保存
        ggsave(
            paste0("output/plot/feature/pseudobulk/", names(subclust_list)[i], "_positive_", gene, ".png"), 
            plot = p_positive, width = 15, height = 10, dpi = 300
        )
    }

    # 各Negative遺伝子のフィーチャープロット（groupで分割）
    for (gene in top_negative_genes) {
        print(paste("Negative Gene:", gene, "feature plot"))
        p_negative <- FeaturePlot(subclust, features = gene, 
                                   cols = c("grey", "blue"), 
                                   ncol = 1, 
                                   split.by = "group") + 
                      theme_classic() + 
                      ggtitle(paste("Negative Gene:", gene))

        # プロットを保存
        ggsave(
            paste0("output/plot/feature/pseudobulk/", names(subclust_list)[i], "_negative_", gene, ".png"), 
            plot = p_negative, width = 15, height = 10, dpi = 300
        )
    }

    # 各Positive遺伝子のバイオリンプロット（groupで分割）
    for (gene in top_positive_genes) {
        print(paste("Positive Gene:", gene, "violin plot"))
        v_positive <- VlnPlot(subclust, features = gene, 
                              split.by = "group",
                              cols = c("grey", "red"), 
                              ncol = 1) + 
                      theme_classic() + 
                      ggtitle(paste("Positive Gene Vln Plot:", gene))

        # プロットを保存
        ggsave(
            paste0("output/plot/violin/pseudobulk/", names(subclust_list)[i], "_positive_", gene, ".png"), 
            plot = v_positive, width = 10, height = 10, dpi = 300
        )
    }

    # 各Negative遺伝子のバイオリンプロット（groupで分割）
    for (gene in top_negative_genes) {
        print(paste("Negative Gene:", gene, "violin plot"))
        v_negative <- VlnPlot(subclust, features = gene, 
                              split.by = "group",
                              cols = c("grey", "blue"), 
                              ncol = 1) + 
                      theme_classic() + 
                      ggtitle(paste("Negative Gene Vln Plot:", gene))

        # プロットを保存
        ggsave(
            paste0("output/plot/violin/pseudobulk/", names(subclust_list)[i], "_negative_", gene, ".png"), 
            plot = v_negative, width = 10, height = 10, dpi = 300
        )
    }
}
