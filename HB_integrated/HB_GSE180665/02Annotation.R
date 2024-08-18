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
library(ComplexHeatmap)

# folder name
folder_name <- "HB_GSE186975"

# harmony_merged_objectを読み込み
seurat_object <- readRDS(paste0("data/", folder_name1, "_", folder_name2, "/harmony_merged_object.rds"))

# 全クラスタの変動遺伝子を抽出
seurat_object_markers <- FindAllMarkers(
    seurat_object,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
)

# marker geneを保存
saveRDS(
    seurat_object_markers,
    file = paste0("data/", folder_name, "/seurat_object_markers.rds")
)


# marker geneを読み込み
# harmony_merged_object_markers <- readRDS(paste0("data/", paste0("Integrated_", folder_name1, "_", folder_name2), "/harmony_merged_object_markers.rds"))

# 各クラスタのDEGの数を確認し分割してリストにまとめる
table(seurat_object_markers$cluster)
deg_list <- split(seurat_object_markers$gene, seurat_object_markers$cluster)

immune_deg_list <- deg_list[c("T_cells", "activated_T_cells", "kupffer", "B_cells", "plasma_cells")]

# リスト内要素の共通項のパターンを集計
common_genes <- make_comb_mat(deg_list)
immune_common_genes <- make_comb_mat(immune_deg_list)

# UpSet Plotを描画
UpSet(
    immune_common_genes,
    set_order = c("T_cells", "activated_T_cells", "kupffer", "B_cells", "plasma_cells"),
)


# 各クラスタの変動遺伝子上位10個を抽出
seurat_object_markers_top10 <- seurat_object_markers %>%
    group_by(cluster) %>%
    top_n(10, avg_log2FC)

seurat_object_markers_top10 %>%
    filter(cluster == "cholangiocyte") %>%
    select(gene)

# DotPlotで可視化
DotPlot(
    seurat_object,
    features = unique(seurat_object_markers_top10$gene),
    group.by = c("seurat_clusters")
) + RotatedAxis() + coord_flip()

# アノテーションのために必要な遺伝子をリストアップ
hepatocyte <- c(
    "AFP", "APOA1", "APOC3", "TTR", "FGB",
    "AHSG", "FABP1", "APOC1", "TF",
    "APOM", "RBP4", "FGA", "FGG", "AMBP",
    "ALB", "VTN", "ORM2", "APOC2", "SLC2A2",
    "HULC", "CYP2D6", "ACOX2", "HFE2"
)
cholangiocyte <- c(
    "KRT8", "KRT18", "TM4SF4", "FXYD2", "C3",
    "SERPING1", "DEFB1", "ANXA4", "KRT7", "ELF3",
    "CYP3A5", "SORBS2", "CLDN4", "CXCL6",
    "CLDN3", "SPP1", "SOX4", "PKHD1", "MTRNR2L8",
    "KRT19", "DSG2", "SOX9", "CFHR2", "CYP2A7"
)
stellate <- c(
    "DCN", "COL3A1", "COL6A2", "COL1A1", "BGN",
    "PTN", "COL1A2", "COLEC11", "COL6A1", "MFAP4",
    "OLFML3", "IGFBP3", "MEST", "COL6A3", "LRRC17"
)
endothelial <- c(
    "RAMP2", "DNASE1L3", "FCN3", "CLDN5", "CAV1",
    "LDB2", "CRHBP", "SDPR", "OIT3", "MYCT1",
    "KDR", "ESAM", "HYAL2", "CLEC1B", "STAB2",
    "AKAP12"
)
T_cells <- c(
    "CD3D", "CD3E", "CD3G", "CD8A", "CD8B",
    "CD4", "TRAC", "TRBC2", "TRGC1",
    "IL2RA", "FOXP3", "CD28", "CTLA4"
)
proliferation <- c(
    "MKI67", "CENPF", "RRM2", "AURKB", "ASPM",
    "ASF1B", "GTSE1", "DLGAP5", "CDCA5", "KIF2C"
)
macrophage <- c(
    "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "LSP1", "CD37",
    "MNDA", "CORO1A", "RETN", "CD48", "RNASE6",
    "HCST", "LYZ"
)
kupffer <- c(
    "CD5L", "TIMD4", "LILRB5", "CETP", "FCGR3A",
    "SLC40A1", "CD163", "MARCO", "FOLR2", "VSIG4",
    "MCOLN1", "LGMN", "MS4A7", "CREG1", "LIPA",
    "C1QC", "C1QB", "C1QA", "CD68"
)
B_cells <- c(
    "CD19", "CD22", "CD79A", "CD79B", "MS4A1",
    "PAX5", "CD40", "BANK1", "TNFRSF13C", "SPIB",
    "VPREB3", "IGHD", "FCRLA", "BCL11A"
)
plasma_cells <- c(
    "IGHG2", "IGHG1", "IGHG3", "TNFRSF17",
    "IGLC2", "IGLL5"
)

features <- c(hepatocyte, cholangiocyte, stellate, endothelial, T_cells, proliferation, kupffer, B_cells, plasma_cells)

# key遺伝子についてDotPlotで可視化
DotPlot(
    seurat_object,
    features = features,
    group.by = "seurat_clusters"
) + RotatedAxis() +
    coord_flip() +
    theme(
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 12),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
    )

ggsave(
    paste0("plot/", folder_name, "/DotPlot_key_genes.png"),
    width = 10, height = 25, dpi = 300, limitsize = FALSE
)



# クラスターに名前をつける
new_cluster_names <- c(
    "13" = "Hepatocyte", "21" = "Hepatoblastoma1", "16" = "Hepatoblastoma2",
    "14" = "Cholangiocyte", "22" = "Stellate", "6" = "Tumor_associated_stellate",
    "9" = "Endothelial", "10" = "Endothelial", "3" = "T_cells", "25" = "T_cells",
    "20" = "proliferation_T_cells", "4" = "Kupffer", "5" = "Kupffer",
    "15" = "B_cells"
)

# クラスター情報を含むメタデータ列の名前を指定します（例: seurat_clusters）
cluster_column <- "seurat_clusters"

# Seuratオブジェクトのクラスターのレベルを取得
cluster_levels <- levels(seurat_object)

# 名前が指定されていないクラスターに "Unknown" を割り当てる
for (level in cluster_levels) {
    if (!(level %in% names(new_cluster_names))) {
        new_cluster_names[level] <- "Unknown"
    }
}

# クラスター名を更新する
seurat_object <- RenameIdents(seurat_object, new_cluster_names)

# 結果を確認
head(Idents(seurat_object))

# meta dataにcell_typeとしてアイデンティティを追加
seurat_object <- AddMetaData(
    seurat_object,
    metadata = Idents(seurat_object),
    col.name = "cell_type"
)

# UMAPをプロット
DimPlot(
    seurat_object,
    reduction = "umap",
    group.by = "cell_type",
    split.by = "tissue2",
    label = TRUE, pt.size = 0.5
)

# feature plot
FeaturePlot(
    seurat_object,
    features = "PAGE4",
    split.by = "tissue2",
    pt.size = 0.5, raster = FALSE
)


# objectにordered_cell_typeを追加
custom_order <- c(
    "hepatocyte", "cholangiocyte", "stellate", "endothelial", "T_cells", "activated_T_cells",
    "kupffer", "B_cells", "plasma_cells"
)

factor(harmony_merged_object$cell_type, levels = custom_order)

harmony_merged_object <- AddMetaData(
    harmony_merged_object,
    metadata = factor(harmony_merged_object$cell_type, levels = custom_order),
    col.name = "ordered_cell_type"
)

### UMAPをプロット
# そのまま表示
DimPlot(
    harmony_merged_object,
    reduction = "umap",
    group.by = "ordered_cell_type",
    label = TRUE, pt.size = 0.5
)
ggsave(
    paste0("plot/Integrated_", folder_name1, "_", folder_name2, "/umap_cell_type.png"),
    width = 10, height = 10, dpi = 300
)

# harmony_merged_objectを保存
saveRDS(
    harmony_merged_object,
    file = paste0("data/", paste0("Integrated_", folder_name1, "_", folder_name2), "/harmony_merged_object.rds")
)
