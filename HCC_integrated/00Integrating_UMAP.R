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
library(harmony)
library(future)
library(sctransform)
library(tidyverse)

# folder name
folder_name <- "HCC_ALL_integrated"
folder_name1 <- "HCC_GSE149614"
folder_name2 <- "HCC_GSE151530"
folder_name3 <- "HCC_GSE124395"
folder_name4 <- "HCC_CELLxGENE_Andrews_2024_Journal_of_Hepatology"

# seurat_objectを読み込み
seurat_object1 <- readRDS(paste0("data/", folder_name1, "/seurat_object.rds"))
seurat_object2 <- readRDS(paste0("data/", folder_name2, "/seurat_object.rds"))
seurat_object3 <- readRDS(paste0("data/", folder_name3, "/seurat_object.rds"))
seurat_object4 <- readRDS(paste0("data/", folder_name4, "/seurat_object.rds"))

# seurat_objectの確認
head(seurat_object1)
seurat_object1
head(seurat_object2)
seurat_object2
head(seurat_object3)
seurat_object3
head(seurat_object4)
seurat_object4

# 各アッセイの遺伝子名を抽出
genes1 <- rownames(seurat_object1[[DefaultAssay(seurat_object1)]])
genes2 <- rownames(seurat_object2[[DefaultAssay(seurat_object2)]])
genes3 <- rownames(seurat_object3[[DefaultAssay(seurat_object3)]])
genes4 <- rownames(seurat_object4[[DefaultAssay(seurat_object4)]])

# 一致しているかどうかを確認
common_genes <- Reduce(intersect, list(genes1, genes2, genes3, genes4))

# 各アッセイの遺伝子名の個数を表示
cat("Genes in seurat_object1:", length(genes1), "\n")
cat("Genes in seurat_object2:", length(genes2), "\n")
cat("Genes in seurat_object3:", length(genes3), "\n")
cat("Genes in seurat_object4:", length(genes4), "\n")
cat("Common genes across all objects:", length(common_genes), "\n")

# Seuratオブジェクトを結合する
# common_genesでフィルタリング
seurat_object1_filtered <- subset(seurat_object1, features = common_genes)
seurat_object2_filtered <- subset(seurat_object2, features = common_genes)
seurat_object3_filtered <- subset(seurat_object3, features = common_genes)
seurat_object4_filtered <- subset(seurat_object4, features = common_genes)

# Seuratオブジェクトをマージ
merged_seurat_object <- merge(
    seurat_object1_filtered,
    y = list(seurat_object2_filtered, seurat_object3_filtered, seurat_object4_filtered),
    add.cell.ids = c(folder_name1, folder_name2, folder_name3, folder_name4),
    project = "HCC_ALL_Integrated"
)

merged_object <- JoinLayers(merged_seurat_object)

# RNA カウントデータを取得
rna_counts <- GetAssayData(merged_object, assay = "RNA", layer = "counts")

# nCount_RNA を再計算 (各細胞の RNA の総カウント数)
nCount_RNA <- colSums(rna_counts)

# nFeature_RNA を再計算 (各細胞の検出された遺伝子数)
nFeature_RNA <- colSums(rna_counts > 0)

# 再計算した値をメタデータに埋め戻す
merged_object@meta.data$nCount_RNA <- nCount_RNA
merged_object@meta.data$nFeature_RNA <- nFeature_RNA

# 計算結果を確認
sum(is.na(merged_object@meta.data$nCount_RNA))
sum(is.na(merged_object@meta.data$nFeature_RNA))

# 新しい列 "split_param" を作成
merged_object$split_param <- NA

# 条件に基づいて "split_param" の値を設定
merged_object$split_param <- case_when(
    merged_object$tissue == "Tumor" ~ "Tumor",
    merged_object$etiology == "Healthy" ~ "Healthy",
    merged_object$etiology %in% c("HBV", "HCV", "HBV,HDV") ~ "Viral",
    merged_object$etiology == "Fatty liver" ~ "Fatty liver",
    merged_object$etiology == "None" ~ "Unknown",
    TRUE ~ NA_character_ # その他の場合はNA
)


# metadataを確認
head(merged_object)
unique(merged_object@meta.data$batch)
unique(merged_object@meta.data$patient)
unique(merged_object@meta.data$tissue)
unique(merged_object@meta.data$stage)
unique(merged_object@meta.data$etiology)
unique(merged_object@meta.data$GSE_ID)
unique(merged_object@meta.data$split_param)

# 確認用
merged_object@meta.data %>%
    filter(split_param == "Tumor") %>%
    select(split_param, etiology, tissue) %>%
    unique()

# merged_object以外のオブジェクトを削除
rm(seurat_object1, seurat_object2, seurat_object3, seurat_object4, seurat_object1_filtered, seurat_object2_filtered, seurat_object3_filtered, seurat_object4_filtered, merged_seurat_object)

# merged_object_RNA_bpという名前でmerged_objectをいったん保存
saveRDS(
    merged_object,
    file = paste0("data/", folder_name, "/merged_object_RNA_bp.rds")
)

# merged_object_RNA_bpを読み込む
merged_object <- readRDS(paste0("data/", folder_name, "/merged_object_RNA_bp.rds"))



### test objectをつかって以下のコードをテスト
## フィルタリング
# ミトコンドリア遺伝子のパーセンテージ計算
merged_object[["percent.mt"]] <- PercentageFeatureSet(merged_object, pattern = "^MT-")

# データのフィルタリング
merged_object <- subset(merged_object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

# 細胞数が少なすぎる(50細胞以下)のbatchを削除
# バッチごとの細胞数を計算
batch_counts <- table(merged_object$batch)

# 細胞数が50を超えるバッチを特定
batches_to_keep <- names(batch_counts[batch_counts > 50])
length(batch_counts)
length(batches_to_keep)

# 選択したバッチのみを含む新しいオブジェクトを作成
merged_object <- subset(merged_object, subset = batch %in% batches_to_keep)
table(merged_object$batch)

# ## テスト用オブジェクトの作成
# set.seed(42)  # 再現性のため

# # 各バッチの細胞数の10%をランダムに選択
# cells_to_keep <- merged_object@meta.data %>%
#   as.data.frame() %>%
#   rownames_to_column("cell_id") %>%
#   group_by(batch) %>%
#   slice_sample(prop = 0.2, replace = FALSE) %>%
#   pull(cell_id)

# # 選択された細胞のみを含む新しいオブジェクトを作成
# test_object <- subset(merged_object, cells = cells_to_keep)

# # バッチごとの細胞数を取得
# batch_counts <- table(test_object$batch)

# # 30細胞以上のバッチを選択
# batches_to_keep <- names(batch_counts[batch_counts > 30])

# # 選択したバッチのみを含む新しいオブジェクトを作成
# test_object_filtered <- subset(test_object, subset = batch %in% batches_to_keep)

# # 結果を確認
# print(table(test_object$batch))
# print(length(test_object$batch))
# print(table(test_object_filtered$batch))
# print(length(test_object_filtered$batch))

# # テストオブジェクトの保存
# saveRDS(test_object_filtered, file = paste0("data/", folder_name, "/test_object.rds"))

# merged_object <- test_object_filtered
# merged_object <- readRDS(paste0("data/", folder_name, "/test_object.rds"))


## 正規化と統合
# layerを分ける
merged_object[["RNA"]] <- split(merged_object[["RNA"]], f = merged_object$batch)
Layers(merged_object)

# SCTransformの実行
merged_object <- SCTransform(merged_object, vars.to.regress = "percent.mt")

# PCAの実行
merged_object <- RunPCA(merged_object)

# UMAPの実行
merged_object <- RunUMAP(merged_object, dims = 1:30)

# 結果を表示
DimPlot(merged_object, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "GSE_ID")

# test folderに保存
ggsave(
    paste0("plot/", folder_name, "/test/test_object_umap.png"),
    width = 20,
    height = 10
)

# harmonyでintegration
merged_object <- IntegrateLayers(
    object = merged_object, 
    method = HarmonyIntegration,
    normalization.method = "SCT"
)
merged_object <- FindNeighbors(merged_object, reduction = "harmony", dims = 1:30)
merged_object <- FindClusters(merged_object, resolution = 0.1)

# UMAPの実行（Harmony reductionまたはPCAを使用）
merged_object <- RunUMAP(merged_object, reduction = "harmony", dims = 1:30)

# UMAPをプロット
DimPlot(merged_object, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "GSE_ID")

# test folderに保存
ggsave(
    paste0("plot/", folder_name, "/test/test_object_umap_harmony.png"),
    width = 20,
    height = 10
)

# 計算後のseurat_objectを保存
saveRDS(
    merged_object,
    file = paste0("data/", folder_name, "/merged_object.rds")
)