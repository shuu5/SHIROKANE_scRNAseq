# libraryの読み込み
library(Seurat)
library(magrittr)
library(glmGamPoi)
library(patchwork)
library(hdf5r)
library(SeuratDisk)
library(SingleR)
library(SingleCellExperiment)
library(scuttle)
library(celldex)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)

# folder name
folder_name <- "HB_GSE180665"

# source scriptの読み込み
source(paste0("src/scripts/load_and_create_seurat.R"))
source(paste0("src/scripts/preprocess1.R"))
source(paste0("src/scripts/SingleRanno.R"))


### データの読み込み、アノテーション
# .h5ad ファイルを .h5seurat 形式に変換
SeuratDisk::Convert(
    source = paste0("data/", folder_name, "/GSE180665.h5ad"),
    dest = paste0("data/", folder_name, "/GSE180665.h5seurat"),
    overwrite = TRUE
)

# .h5seurat ファイルを読み込んで Seurat オブジェクトを作成
sc <- SeuratDisk::LoadH5Seurat(file = paste0("data/", folder_name, "/GSE180665.h5seurat"))
sc
head(sc)
unique(sc@meta.data$sample_group)

metadata <- sc@meta.data

head(metadata)

# sample_groupがpdxもしくはtumorのデータをtissue列に"tumor"として追加
metadata$tissue <- ifelse(metadata$sample_group %in% c("pdx", "tumor"), "tumor", "normal")

# Seuratオブジェクトにメタデータを追加
sc <- AddMetaData(sc, metadata)

# UMAPを描画
DimPlot(sc, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "tissue")

# feature plot
FeaturePlot(sc, features = "PAGE4", split.by = "tissue", pt.size = 0.5, raster = FALSE)
VlnPlot(sc, features = "PAGE4", split.by = "tissue")



hepatocyte <- c(
    "AFP", "APOA1", "APOC3", "TTR", "FGB",
    "AHSG", "FABP1", "APOC1", "TF",
    "APOM", "RBP4", "FGA", "FGG", "AMBP",
    "ALB", "VTN", "ORM2", "APOC2", "SLC2A2",
    "HULC", "CYP2D6", "ACOX2", "HFE2"
)

T_cells <- c(
    "CD3D", "CD3E", "CD3G", "CD8A", "CD8B",
    "CD4", "TRAC", "TRBC2", "TRGC1",
    "IL2RA", "FOXP3", "CD28", "CTLA4"
)

FeaturePlot(sc, features = "CD3D", split.by = "tissue", pt.size = 0.5, raster = FALSE)
