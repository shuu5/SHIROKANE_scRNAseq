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
library(metap)


# folder name
folder_name <- "HCC_ALL_integrated"

# merged_objectを読み込み
merged_object <- readRDS(paste0("data/", folder_name, "/merged_object.rds"))


# 全クラスタの変動遺伝子を抽出
merged_object <- PrepSCTFindMarkers(merged_object)
merged_object_markers <- FindAllMarkers(
    merged_object,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
)

# adj p valが0のものを除外
merged_object_markers <- merged_object_markers %>%
    filter(p_val_adj > 0)

# marker geneを保存
saveRDS(
    merged_object_markers,
    file = paste0("data/", folder_name, "/merged_object_markers.rds")
)

# marker geneを読み込み
# merged_object_markers <- readRDS(paste0("data/", folder_name, "/merged_object_markers.rds"))


# 各クラスタの変動遺伝子上位10個を抽出
merged_object_markers_top10 <- merged_object_markers %>%
    group_by(cluster) %>%
    top_n(10, avg_log2FC)

print(merged_object_markers_top10, n=150)

merged_object_markers_top10 %>%
    filter(cluster == "8")

# DotPlotで可視化
DotPlot(
    merged_object,
    features = unique(merged_object_markers_top10$gene),
    group.by = c("seurat_clusters")
) + RotatedAxis() + coord_flip() + theme_classic()

# 画像を保存
ggsave(
    paste0("plot/", folder_name, "/Marker_top10_DotPlot_UMAP.png"),
    width = 10, height = 20, dpi = 300
)


# 各クラスタの変動遺伝子からクラスターに名前を付ける
new_cluster_names <- c(
    "1" = "Hepatocyte_Mitochondrial", 
    "2" = "Kupffer_Cell", 
    "3" = "T_cell_1",
    "4" = "Hepatocyte_Metabolic", 
    "5" = "T_cell_2", 
    "6" = "Endothelial",
    "7" = "Cholangiocyte", 
    "8" = "Hepatic_Stellate_Cell",
    "9" = "Stressed_Hepatocyte",
    "10" = "Liver_Resident_NK_Cell", 
    "11" = "Liver_Progenitor_Cell",
    "12" = "Activated_T_Cell",
    "13" = "Sinusoidal_Endothelial_Cell", 
    "14" = "Liver_Specific_Sinusoidal_Endothelial_Cell",
    "15" = "Drug_Metabolizing_Hepatocyte"
)
names(new_cluster_names) <- levels(merged_object)

# アイデンティティを更新する
merged_object <- RenameIdents(merged_object, new_cluster_names)

# meta dataにcell_typeとしてアイデンティティを追加
merged_object <- AddMetaData(
    merged_object,
    metadata = Idents(merged_object),
    col.name = "cell_type"
)

### UMAPをプロット
DimPlot(
    merged_object,
    reduction = "umap",
    group.by = "cell_type",
    label = TRUE, pt.size = 0.5
)

ggsave(
    paste0("plot/", folder_name, "/umap_cell_type.png"),
    width = 15, height = 10, dpi = 300
)

# merged_objectを保存
saveRDS(
    merged_object,
    file = paste0("data/", folder_name, "/merged_object.rds")
)




# DotPlotで可視化
DotPlot(
    merged_object,
    features = hepatocyte,
    group.by = c("cell_type")
) + RotatedAxis() + theme_classic()

DotPlot(
    merged_object,
    features = cholangiocyte,
    group.by = c("cell_type")
) + RotatedAxis() + theme_classic()

DotPlot(
    merged_object,
    features = stellate,
    group.by = c("cell_type")
) + RotatedAxis() + theme_classic()

DotPlot(
    merged_object,
    features = endothelial,
    group.by = c("cell_type")
) + RotatedAxis() + theme_classic()

DotPlot(
    merged_object,
    features = T_cells,
    group.by = c("cell_type")
) + RotatedAxis() + theme_classic()

DotPlot(
    merged_object,
    features = proliferation,
    group.by = c("cell_type")
) + RotatedAxis() + theme_classic()

DotPlot(
    merged_object,
    features = macrophage,
    group.by = c("cell_type")
) + RotatedAxis() + theme_classic()

DotPlot(
    merged_object,
    features = kupffer,
    group.by = c("cell_type")
) + RotatedAxis() + theme_classic()

DotPlot(
    merged_object,
    features = B_cells,
    group.by = c("cell_type")
) + RotatedAxis() + theme_classic()

DotPlot(
    merged_object,
    features = plasma_cells,
    group.by = c("cell_type")
) + RotatedAxis() + theme_classic()




### 既知のmarker gene setを利用する方法
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
    "CD24", "CYP3A5", "SORBS2", "CLDN4", "CXCL6",
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
    "CD4", "TRAC", "TRBC1", "TRBC2", "TRGC1",
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
    "JCHAIN", "IGHG2", "IGHG1", "IGHG3", "TNFRSF17",
    "IGLC2", "IGHGP", "IGLL5"
)

# features <- c(hepatocyte, cholangiocyte, stellate, endothelial, T_cells, proliferation, macrophage, kupffer, B_cells, plasma_cells)

# # key遺伝子についてDotPlotで可視化
# DotPlot(
#     merged_object,
#     features = features,
#     group.by = "seurat_clusters"
# ) + RotatedAxis() +
#     coord_flip() +
#     theme(
#         axis.text.x = element_text(size = 18),
#         axis.text.y = element_text(size = 12),
#         panel.background = element_rect(fill = "white"),
#         plot.background = element_rect(fill = "white", colour = "white"),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank()
#     )

# ggsave(
#     paste0("plot/", folder_name, "/DotPlot_key_genes.png"),
#     width = 10, height = 25, dpi = 300, limitsize = FALSE
# )

# DimPlot(merged_object, reduction = "umap", split.by = "tissue", label = TRUE, pt.size = 0.5)

# # クラスターに名前をつける
# new_cluster_names <- c(
#     "1" = "T_cells", "2" = "hepatocyte", "3" = "kupffer",
#     "4" = "cholangiocyte", "5" = "endothelial", "6" = "Activated_T_cells",
#     "7" = "stellate", "8" = "plasma_cells", "9" = "B_cells"
# )
# names(new_cluster_names) <- levels(merged_object)

# # アイデンティティを更新する
# merged_object <- RenameIdents(merged_object, new_cluster_names)

# # meta dataにcell_typeとしてアイデンティティを追加
# merged_object <- AddMetaData(
#     merged_object,
#     metadata = Idents(merged_object),
#     col.name = "cell_type"
# )

# # objectにordered_cell_typeを追加
# custom_order <- c(
#     "hepatocyte", "cholangiocyte", "stellate", "endothelial", "T_cells",
#     "kupffer", "B_cells", "plasma_cells", "Activated_T_cells"
# )

# factor(merged_object$cell_type, levels = custom_order)

# merged_object <- AddMetaData(
#     merged_object,
#     metadata = factor(merged_object$cell_type, levels = custom_order),
#     col.name = "ordered_cell_type"
# )

# ### UMAPをプロット
# # そのまま表示
# DimPlot(
#     merged_object,
#     reduction = "umap",
#     group.by = "ordered_cell_type",
#     label = TRUE, pt.size = 0.5
# )

# ggsave(
#     paste0("plot/", folder_name, "/umap_ordered_cell_type.png"),
#     width = 10, height = 10, dpi = 300
# )

# # merged_objectを保存
# saveRDS(
#     merged_object,
#     file = paste0("data/", folder_name, "/merged_object.rds")
# )