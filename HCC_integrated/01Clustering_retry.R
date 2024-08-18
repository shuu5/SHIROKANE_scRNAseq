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

# folder name
folder_name <- "HCC_ALL_integrated"

# merged_objectを読み込み
merged_object <- readRDS(paste0("data/", folder_name, "/merged_object.rds"))

# 条件を変えて再クラスタリング
res <- 0.1
merged_object <- FindClusters(
    merged_object,
    resolution = res,
    algorithm = 4,
    method = "igraph"
)

# 可視化
DimPlot(merged_object, reduction = "umap", label = TRUE, pt.size = 0.5)
head(merged_object)

# plotを保存
ggsave(
    paste0("plot/", folder_name, "/umap_res", res, ".png"),
    width = 10, height = 10, dpi = 300
)

# 計算後のseurat_objectを保存
saveRDS(
    merged_object,
    file = paste0("data/", folder_name, "/merged_object.rds")
)
