library(Seurat)
library(tidyverse)

# data/processed/subclusters/内のsubclusteringしたオブジェクトを読み込み
subclust_list <- lapply(list.files("data/processed/subclusters/", pattern = "*.rds", full.names = TRUE), readRDS)
names(subclust_list) <- gsub("data/processed/subclusters/(.*)_subclust.rds", "\\1", list.files("data/processed/subclusters/", pattern = "*.rds"))
names(subclust_list) <- gsub("_subclust.rds", "", names(subclust_list))

# 特に見たいgenesの組み合わせのリスト
tar_genes_list <- list(
  c("ARID5A", "KRT13"), 
  c("ARID5A", "TGM2"), 
  c("ARID5A", "IL6"),
  c("ARID5A", "ARID5B"),
  c("ARID5B", "KRT13"), 
  c("ARID5B", "TGM2"), 
  c("ARID5B", "IL6"),
  c("GAPDH", "MKI67")
  )
total_genes <- unique(unlist(tar_genes_list))



### 各subclust_listのオブジェクトに対して特定の2遺伝子のFeatureScatterをplotする
lapply(tar_genes_list, function(tar_genes) {
  lapply(names(subclust_list), function(cell_type) {
    print(paste("plotting:", cell_type, tar_genes[1], tar_genes[2]))

    subset_obj <- subclust_list[[cell_type]]

    # tar_genesがsubset_objに存在するか確認
    if (!all(tar_genes %in% rownames(subset_obj))) {
      warning(paste("Skipping:", cell_type, "because one or both genes are missing:", tar_genes[1], tar_genes[2]))
      return(NULL)  # 次に進む
    }

    FeatureScatter(subset_obj, feature1 = tar_genes[1], feature2 = tar_genes[2]) +
    theme_classic()

    ggsave(paste0("output/plot/scatter/subclusters/", cell_type, "_", tar_genes[1], "_", tar_genes[2], ".png"), width = 10, height = 10)
  })
})


### RNAのcountsを使ってFeatureScatterをplotする
# すべてのオブジェクトのDefaultAssayを"RNA"に設定
for(i in 1:length(subclust_list)) {
  DefaultAssay(subclust_list[[i]]) <- "RNA"
}

lapply(tar_genes_list, function(tar_genes) {
  lapply(names(subclust_list), function(cell_type) {
    print(paste("plotting:", cell_type, tar_genes[1], tar_genes[2]))

    subset_obj <- subclust_list[[cell_type]]

    # tar_genesがsubset_objに存在するか確認
    if (!all(tar_genes %in% rownames(subset_obj))) {
      warning(paste("Skipping:", cell_type, "because one or both genes are missing:", tar_genes[1], tar_genes[2]))
      return(NULL)  # 次に進む
    }

    FeatureScatter(subset_obj, feature1 = tar_genes[1], feature2 = tar_genes[2], slot = "counts") +
    theme_classic()

    ggsave(paste0("output/plot/scatter/subclusters/", cell_type, "_", tar_genes[1], "_", tar_genes[2], "_RNA.png"), width = 10, height = 10)
  })
})