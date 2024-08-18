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
library(dplyr)

# project名
project_name <- "GC_PRJCA004201"

### データの読み込み、アノテーション
# ファイルのパスを指定
matrix_file <- paste0("data/", project_name, "/matrix.mtx.gz")
barcodes_file <- paste0("data/", project_name, "/barcodes.tsv.gz")
features_file <- paste0("data/", project_name, "/features.tsv.gz")

# ファイルを読み込み、Seuratオブジェクトを作成
seurat_object <- CreateSeuratObject(counts = Read10X(data.dir = paste0("data/", project_name)), project = project_name)

# データの前処理
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

# メタデータファイルを読み込む
metadata <- read.table(paste0("data/", project_name, "/OMIX001073-20-04.tsv"), header = TRUE, sep = "\t", row.names = 1, fill = TRUE)

# metadataからbatch, patient, tissueを抽出
metadata_selected <- metadata[, c("batch", "patient", "tissue")]

# seurat_objectにmetadata_selectedを追加
seurat_object <- AddMetaData(seurat_object, metadata_selected)

head(seurat_object)

# もう一つのmetadataを読み込む
metadata_2 <- read.csv("data/GC_PRJCA004201/metadata.csv", header = TRUE)
metadata_2

# 必要な情報を抜き出す
metadata_2_selected <- metadata_2[, c("Diagnosis", "Lauren", "H.pylori", "Cancer_location", "Patient.ID")]

# Patient.IDをキーとしてmetadata_2_selectedをseurat_objectに追加
rownames(metadata_2_selected) <- metadata_2_selected$Patient.ID
metadata_2_selected <- metadata_2_selected[, -which(names(metadata_2_selected) == "Patient.ID")]

# seurat_objectのpatient列をキーにしてmetadata_2_selectedを追加
seurat_object <- AddMetaData(seurat_object, metadata_2_selected[seurat_object$patient, ])

# 結果を確認
head(seurat_object)

# 結果を保存
saveRDS(seurat_object, file = "data/GC_PRJCA004201/seurat_object.rds")
