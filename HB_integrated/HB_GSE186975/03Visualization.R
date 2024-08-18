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
folder_name <- "HB_GSE186975"

# source scriptの読み込み
source(paste0("src/", folder_name, "/scripts/visualizeTargetGene.R"))

# seurat_objectを読み込み
seurat_object <- readRDS(paste0("data/", folder_name, "/seurat_object_UMAP.rds"))


# メタデータを取得
metadata <- seurat_object@meta.data

# 新しいtissue2列を作成
metadata$tissue2 <- ifelse(metadata$tissue %in% c("H", "N"), "N", "T")

# メタデータをSeuratオブジェクトに戻す
seurat_object@meta.data <- metadata

# 結果を確認
head(seurat_object@meta.data)




# 可視化
DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5)
ggsave(paste0("plot/", folder_name, "/umap_plot.png"), width = 10, height = 10, dpi = 300)

# tissueでsplit
DimPlot(seurat_object, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "tissue2")
ggsave(paste0("plot/", folder_name, "/umap_plot_cell_type.png"), width = 14, height = 10, dpi = 300)



# 特定の遺伝子の発現量を可視化
target_gene <- "YAP1"
VlnPlot(seurat_object, features = target_gene)
FeaturePlot(seurat_object, features = target_gene, pt.size = 0.5, raster = FALSE)
