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
project_name <- "GC_GSE150290"

# ディレクトリのリストを取得
base_dir <- file.path("data", project_name)
dirs <- list.dirs(path = base_dir, full.names = TRUE, recursive = FALSE)

# 各ディレクトリ内の10xデータを読み込む
seurat_list <- list()
for (dir in dirs) {
    # ディレクトリ名を取得
    dir_name <- basename(dir)

    # 10xデータのファイルパスを指定
    data_dir <- file.path(dir)

    # 10xデータを読み込む
    seurat_obj <- Read10X(data.dir = data_dir)

    # Seuratオブジェクトを作成
    seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = dir_name)

    # ミトコンドリア遺伝子のパーセンテージ計算
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

    # データのフィルタリング
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

    # 必要に応じて、Seuratオブジェクトをリストに保存するなどの処理を追加
    seurat_list[[dir_name]] <- seurat_obj

    # デバッグ用にメッセージを表示
    message("Loaded data for: ", dir_name)
    message("Number of cells: ", ncol(seurat_obj))
    message("Number of genes: ", nrow(seurat_obj))
}

# seurat_listに保存された各seurat_objのcell_nameを変更
for (dir_name in names(seurat_list)) {
    seurat_obj <- seurat_list[[dir_name]]

    # セル名を変更
    cell_names <- colnames(seurat_obj)
    new_cell_names <- paste0(dir_name, "_", cell_names)
    colnames(seurat_obj) <- new_cell_names

    # 変更をリストに反映
    seurat_list[[dir_name]] <- seurat_obj
}

# 結果の確認
seurat_list
head(seurat_list[[11]])

# すべてのSeuratオブジェクトをマージ
merged_seurat <- Reduce(function(x, y) merge(x, y), seurat_list)
print(merged_seurat)

# layerを統合
merged_seurat <- JoinLayers(merged_seurat)
print(merged_seurat)
head(merged_seurat)





# metadataの読み込み
metadata <- read.csv(file.path(base_dir, "metadata.csv"), row.names = 1)
metadata

# metadataのうち必要な列を抜き出す
metadata_df <- metadata[, c("Diagnosis", "Lauren", "Atrophy", "H.pylori", "Cancer_location", "MSI", "EBV")]
metadata_df

# Patient1からPatient24までのデータを2行に複製し、名前を変更
metadata_df_expanded <- metadata_df[rep(1:24, each = 2), ]
rownames(metadata_df_expanded) <- paste0("Pat", sprintf("%02d", rep(1:24, each = 2)), "-", rep(c("A", "B"), 24))

# Pat21-Bを削除
metadata_df_expanded <- metadata_df_expanded[rownames(metadata_df_expanded) != "Pat21-B", ]

# Patient25以降のデータを追加し、名前を変更
additional_metadata <- metadata_df[25:nrow(metadata_df), ]
rownames(additional_metadata) <- paste0("Pat", sprintf("%02d", 25:nrow(metadata_df)), "-A")
metadata_df_expanded <- rbind(metadata_df_expanded, additional_metadata)

# 新しい列を追加
metadata_df_expanded$batch <- rownames(metadata_df_expanded)
metadata_df_expanded$patient <- paste0("patient", as.numeric(gsub("Pat|-[AB]", "", rownames(metadata_df_expanded))))
metadata_df_expanded$tissue <- ifelse(grepl("-A$", rownames(metadata_df_expanded)), "Normal", "Tumor")

# 結果の確認
metadata_df_expanded


# orig.ident列を取得
orig_ident <- merged_seurat@meta.data$orig.ident

# metadata_df_expandedの行名と一致するようにメタデータを追加
metadata_to_add <- metadata_df_expanded[orig_ident, ]
head(metadata_to_add)

# メタデータを追加
merged_seurat <- AddMetaData(object = merged_seurat, metadata = metadata_to_add)

# 結果の確認
head(merged_seurat@meta.data)

# 結果の保存
saveRDS(merged_seurat, file = file.path(base_dir, "merged_seurat.rds"))
