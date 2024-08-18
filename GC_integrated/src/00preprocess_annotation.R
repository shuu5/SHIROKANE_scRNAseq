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
project_name <- "GC_integrated"
project_name_1 <- "GC_PRJCA004201"
project_name_2 <- "GC_GSE150290"

# それぞれのデータを読み込み
seurat_object_1 <- readRDS(file.path("data", project_name_1, "seurat_object.rds"))
seurat_object_2 <- readRDS(file.path("data", project_name_2, "merged_seurat.rds"))

# seurat_object_1のtissueがBloodの細胞を取り除く
seurat_object_1 <- subset(seurat_object_1, subset = tissue != "Blood")

# tissue列の修正
seurat_object_2@meta.data <- seurat_object_2@meta.data %>%
    mutate(tissue = ifelse(patient %in% paste0("patient", 1:24) & tissue == "Normal", "Paratumor", tissue))

# 確認
seurat_object_1
seurat_object_2
head(seurat_object_1)
head(seurat_object_2)
unique(seurat_object_1@meta.data$tissue)
unique(seurat_object_2@meta.data$tissue)
table(seurat_object_1@meta.data$tissue)
table(seurat_object_2@meta.data$tissue)

# Seuratオブジェクトをマージ
merged_seurat_object <- merge(
    seurat_object_1,
    y = seurat_object_2,
    add.cell.ids = c(project_name_1, project_name_2),
    project = project_name
)

# 結果の確認
merged_object
head(merged_object)
table(merged_object@meta.data$H.pylori)

# 保存
saveRDS(merged_object, file.path("data", project_name, "merged_seurat.rds"))