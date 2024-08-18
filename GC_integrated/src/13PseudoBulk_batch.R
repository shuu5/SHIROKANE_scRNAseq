library(Seurat)
library(tidyverse)
library(presto)
library(ggrepel)

# data/processed/subclusters/内のsubclusteringしたオブジェクトを読み込み
subclust_list <- lapply(list.files("data/processed/subclusters/", pattern = "*.rds", full.names = TRUE), readRDS)
names(subclust_list) <- gsub("data/processed/subclusters/(.*)_subclust.rds", "\\1", list.files("data/processed/subclusters/", pattern = "*.rds"))
names(subclust_list) <- gsub("_subclust.rds", "", names(subclust_list))

# 特に見たいgenes
tar_genes <- c("ARID5A", "ARID5B", "KRT13", "TGM2", "AXL", "MKI67", "IL6")
vs_genes_list <- list(
  c("ARID5A", "ARID5B"),
  c("ARID5A", "KRT13"),
  c("ARID5A", "TGM2"),
  c("ARID5A", "IL6")
)

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

### subclust_listの各オブジェクトのbatchごとの平均発現量を計算しplotする
for (i in seq_along(subclust_list)) {
    subclust <- subclust_list[[i]]
    
    # 各サブクラスタの平均発現値を計算
    agg_subclust <- AggregateExpression(subclust, group.by = c("batch", "group", "Lauren"), return.seurat = TRUE)

    # agg_subclustの確認
    agg_subclust
    
    # データフレームに変換
    expression_data <- as.data.frame(agg_subclust@assays$SCT$data)
    expression_data$gene <- rownames(expression_data)
    expression_data_long <- pivot_longer(expression_data, -gene, names_to = "sample", values_to = "expression")

    # sampleからbatchとgroupとLaurenを抽出
    expression_data_long <- expression_data_long %>%
        mutate(batch = sub("_.*", "", sample), 
               group = sub(".*_(.*)_.*", "\\1", sample), 
               Lauren = sub(".*_(.*)", "\\1", sample))

    # tar_genesの中でsubclustに含まれている遺伝子をpresent_genesに格納
    present_genes <- intersect(tar_genes, rownames(subclust))

    # 各groupにおける特定のgeneの発現を比較するボックスプロットを作成
    for (tar_gene in present_genes) {
      print(paste0(tar_gene, "_boxplot.png processing..."))
      ggplot(expression_data_long %>% filter(gene == tar_gene), aes(x = group, y = expression)) + 
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(width = 0.2, size = 1) + 
          theme_classic() +
          ggtitle(paste("Gene:", tar_gene)) +
          labs(fill = "Batch")  # バッチを色分け
      ggsave(paste0("output/plot/boxplot/pseudobulk/", names(subclust_list)[i], "_", tar_gene, "_boxplot.png"), width = 6, height = 8)
    }

    # 各group毎に二つの遺伝子を比較するscatter plotを作成
    for (vs_genes in vs_genes_list) {
      if (all(vs_genes %in% present_genes)) {
        print(paste0(vs_genes[1], "_vs_", vs_genes[2], "_scatterplot.png processing..."))

        gene1 = expression_data_long %>% filter(gene == vs_genes[1]) %>% pull(expression)
        gene2 = expression_data_long %>% filter(gene == vs_genes[2]) %>% pull(expression)
        lauren = expression_data_long %>% filter(gene == vs_genes[1]) %>% pull(Lauren)

        ggplot(NULL, aes(x = gene1, y = gene2, color = lauren)) +
            geom_point() +
            theme_classic() +
            labs(x = vs_genes[1], y = vs_genes[2], color = "Lauren")
        ggsave(paste0("output/plot/scatter/pseudobulk/", names(subclust_list)[i], "_", vs_genes[1], "_vs_", vs_genes[2], "_scatterplot.png"), width = 10, height = 10)
      } else {
        print(paste("Skipping scatter plot for", vs_genes[1], "and", vs_genes[2], "as one or both genes are not present."))
      }
    }
}