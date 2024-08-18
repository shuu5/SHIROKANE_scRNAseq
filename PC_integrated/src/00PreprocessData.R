# Seuratパッケージの読み込み
library(Seurat)
library(tidyverse)

# RDSファイルの読み込み
counts_data <- readRDS("data/raw/GSE210347_counts.rds")

# metaデータのよみこみ(GSE210347_meta.txt, タブ区切り)
meta_data <- read.table("data/raw/GSE210347_meta.txt", header = TRUE, sep = "\t")

# meta_data_pdacを取得
meta_data_pdac <- meta_data %>%
    filter(tissue == "PDAC")

# batch列の作成
meta_data_pdac$batch <- ifelse(grepl("PUMCH", meta_data_pdac$SampleID),
                               sub(".*_(PUMCH_\\w+).*", "\\1", meta_data_pdac$SampleID),
                               "eliminated")

# eliminatedの削除
meta_data_pdac <- meta_data_pdac %>%
    filter(batch != "eliminated")

# pdacのcellnameを取得
pdac_cellname <- meta_data_pdac$cellname

# counts_dataからPDACの細胞を抜き出す
pdac_counts_data <- counts_data[, pdac_cellname]

# seurat_obj_pdacを作成
seurat_obj_pdac <- CreateSeuratObject(counts = pdac_counts_data, project = "GSE210347_Project")

# meta_data_pdacをseurat_obj_pdacに追加
seurat_obj_pdac@meta.data <- meta_data_pdac

# 臨床データの読み込み(data/raw/GSE210437_PC_patient_data.csv)
clinical_data_pdac <- read.csv("data/raw/GSE210437_PC_patient_data.csv", skip = 1) 

# Number列の値を変更
clinical_data_pdac$Number <- gsub("-", "_", clinical_data_pdac$Number)

# seurat_obj_pdacのメタデータとclinical_data_pdacをマージ
merged_meta_data <- merge(seurat_obj_pdac@meta.data, clinical_data_pdac, by.x = "batch", by.y = "Number", all.x = TRUE)

# マージしたメタデータをseurat_obj_pdacに追加
seurat_obj_pdac@meta.data <- merged_meta_data

# カウントマトリックスの列名を確認
count_colnames <- colnames(seurat_obj_pdac[["RNA"]]$counts)
head(count_colnames)

# メタデータの行名を確認
meta_rownames <- rownames(seurat_obj_pdac@meta.data)
head(meta_rownames)

# メタデータの行名をカウントマトリックスの列名に合わせる
rownames(seurat_obj_pdac@meta.data) <- count_colnames

# 不要なオブジェクトの削除
rm(counts_data, meta_data, meta_data_pdac, pdac_cellname, pdac_counts_data, clinical_data_pdac, merged_meta_data)

# メモリの解放
gc()

# 結果の確認
print(seurat_obj_pdac)
print(head(seurat_obj_pdac@meta.data))
print(table(seurat_obj_pdac@meta.data$group, seurat_obj_pdac@meta.data$celltype))
print(unique(seurat_obj_pdac@meta.data$SampleID))
print(unique(seurat_obj_pdac@meta.data$batch))
print(table(seurat_obj_pdac$batch))

# seurat_obj_pdacを保存
saveRDS(seurat_obj_pdac, "data/processed/seurat_obj_pdac.rds")
