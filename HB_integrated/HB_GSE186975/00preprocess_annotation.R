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
source(paste0("src/scripts/load_and_create_seurat.R"))
source(paste0("src/scripts/preprocess1.R"))
source(paste0("src/scripts/SingleRanno.R"))


### データの読み込み、アノテーション
# folder内のすべてのファイル名を取得する
file_names <- list.files(paste0("data/", folder_name), pattern = "dge.txt", full.names = TRUE)
file_names

# ファイル名からメタデータを取得
# "GSMxxxxxxx"と最後の英字を抽出する関数
extract_gsm_and_suffix <- function(file_name) {
    gsm <- stringr::str_extract(file_name, "GSM\\d{7}")
    suffix <- stringr::str_extract(file_name, "(?<=DB_HB\\d{3})[A-Z]")
    return(c(gsm, suffix))
}

# 各ファイル名から抽出
results <- t(sapply(file_names, extract_gsm_and_suffix))

# 結果をデータフレームとして表示
results_df <- as.data.frame(results, stringsAsFactors = FALSE)
colnames(results_df) <- c("GSM", "Suffix")
print(results_df)


# ファイルを読み込み、Seuratオブジェクトを作成
seurat_object <- load_and_create_seurat_object(paste0("data/", folder_name, "/GSM5664587_DB_HB001H_Pool_A_B_dge.txt"), sep = "\t")
seurat_object
head(seurat_object)

# Seuratオブジェクトにメタデータを追加
seurat_object@meta.data$batch <- results_df$GSM[1]
seurat_object@meta.data$tissue <- results_df$Suffix[1]
head(seurat_object)



library(Seurat)

# メタデータの抽出関数（以前のコードと同じ）
extract_gsm_and_suffix <- function(file_name) {
    gsm <- stringr::str_extract(file_name, "GSM\\d{7}")
    suffix <- stringr::str_extract(file_name, "(?<=DB_HB\\d{3})[A-Z]")
    return(c(gsm, suffix))
}

# メタデータの抽出
results <- t(sapply(file_names, extract_gsm_and_suffix))
results_df <- as.data.frame(results, stringsAsFactors = FALSE)
colnames(results_df) <- c("GSM", "Suffix")

# Seuratオブジェクトの作成とメタデータの追加
seurat_list <- list()
for (i in 1:length(file_names)) {
    seurat_object <- load_and_create_seurat_object(file_names[i], sep = "\t")
    seurat_object@meta.data$batch <- results_df$GSM[i]
    seurat_object@meta.data$tissue <- results_df$Suffix[i]

    # セル名を一意にする
    cell_names <- colnames(seurat_object)
    new_cell_names <- paste0(results_df$GSM[i], "_", cell_names)
    colnames(seurat_object) <- new_cell_names

    seurat_list[[i]] <- seurat_object
}

# 結果の表示
head(seurat_list[[6]]@meta.data)

# Merge all Seurat objects into one
merged_seurat_object <- merge(x = seurat_list[[1]], y = seurat_list[-1])

# Ensure that all counts are in the "RNA" assay
DefaultAssay(merged_seurat_object) <- "RNA"

# Remove unnecessary layers (assays)
Assays(merged_seurat_object) <- "RNA"


# JoinLayers関数を使用してレイヤーを一つに統合
merged_seurat_object <- JoinLayers(merged_seurat_object)

# Check the result
print(merged_seurat_object)
head(merged_seurat_object@meta.data)

# Seuratオブジェクトを保存
saveRDS(merged_seurat_object, paste0("data/", folder_name, "/merged_seurat_object.rds"))


# データの前処理
seurat_object <- process_seurat_data(merged_seurat_object)


# 予め用意された参照データセットをロード
ref <- HumanPrimaryCellAtlasData()
class(ref)
# SingleRを用いてアノテーションを行う
seurat_object <- annotate_seurat_with_singleR(seurat_object)

# 別のリファレンスでアノテーションを行う(手動でダウンロードしてSingleCellExperiment形式に変換する)

# オブジェクトの検証
head(seurat_object)
head(seurat_object@meta.data)
ncol(seurat_object)
nrow(seurat_object)
ncol(seurat_object@meta.data)
nrow(seurat_object@meta.data)
validObject(seurat_object)

# seurat_objectを一度保存
saveRDS(seurat_object, paste0("data/", folder_name, "/seurat_object.rds"))
