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
folder_name <- "HCC_ALL_integrated"

# harmony_merged_objectを読み込み
merged_object <- readRDS(paste0("data/", folder_name, "/merged_object.rds"))

# 各meta.dataの情報を確認
merged_object
head(merged_object)
table(merged_object@meta.data$tissue)
table(merged_object@meta.data$split_param)
table(merged_object@meta.data$etiology)
table(merged_object@meta.data$cell_type)
table(merged_object@meta.data$ordered_cell_type)

# tissueでわけたUMAPを可視化
DimPlot(
    merged_object,
    reduction = "umap",
    group.by = "cell_type",
    split.by = "tissue",
    label = TRUE, pt.size = 0.5
)

ggsave(
    paste0("plot/", folder_name, "/umap_celltype_split_by_tissue.png"),
    width = 20, height = 10, dpi = 300
)

# split_paramでわけたUMAPを可視化
DimPlot(
    merged_object,
    reduction = "umap",
    group.by = "cell_type",
    split.by = "split_param",
    label = TRUE, pt.size = 0.5
)

ggsave(
    paste0("plot/", folder_name, "/umap_celltype_split_param.png"),
    width = 30, height = 10, dpi = 300
)

# etiologyでわけたUMAPを可視化
DimPlot(
    merged_object,
    reduction = "umap",
    group.by = "cell_type",
    split.by = "etiology",
    label = TRUE, pt.size = 0.5
)

ggsave(
    paste0("plot/", folder_name, "/umap_celltype_etiology.png"),
    width = 40, height = 10, dpi = 300
)


### target geneに対するplot
# target geneをまとめる
total_genes <- c(
    "KCNN4", "ATG9B", "STAC3", "HOXA3", "HOXA2", 
    "ACAP1", "TMC8", "CLEC2D", "TNFRSF25", "IFITM1", 
    "IFITM2", "IFITM3", "CD74", "CD177", "CAR1", "FABP2"
    )

# DotPlot
DotPlot(
    merged_object,
    features = total_genes,
    group.by = "cell_type"
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
    paste0("plot/", folder_name, "/DotPlot_target_genes.png"),
    width = 12, height = 8, dpi = 300, limitsize = FALSE
)

# DotPlotをsplit.by = "tissue"で分ける
DotPlot(
    merged_object,
    features = total_genes,
    group.by = "cell_type",
    split.by = "tissue",
    cols = c("blue", "blue")  # 色を指定
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
    paste0("plot/", folder_name, "/DotPlot_target_genes_tissue.png"),
    width = 20, height = 8, dpi = 300, limitsize = FALSE
)


# 各種plotをtotal_genesに対して施行
for (i in 1:length(total_genes)) {
    # そのままのfeatureplotとviolinplot
    Feature <- FeaturePlot(
        merged_object,
        features = total_genes[i],
        pt.size = 0.5, label = TRUE, order = TRUE
    )
    Vln <- VlnPlot(
        merged_object,
        features = total_genes[i],
        pt.size = 0
    )
    Feature + Vln
    ggsave(
        paste0("plot/", folder_name, "/total_target/Feature_Vln_", total_genes[i], ".png"),
        width = 10, height = 20, dpi = 300
    )

    # split.by = "tissue"でNormalとTumorで分ける
    Feature_tissue <- FeaturePlot(
        merged_object,
        features = total_genes[i],
        pt.size = 0.5, label = TRUE, order = TRUE,
        split.by = "tissue"
    )

    Vln_tissue <- VlnPlot(
        merged_object,
        features = total_genes[i],
        pt.size = 0,
        split.by = "tissue"
    )

    plot_grid(Feature_tissue, Vln_tissue, ncol = 1, nrow = 2)
    ggsave(
        paste0("plot/", folder_name, "/total_target/Feature_Vln_tissue_", total_genes[i], ".png"),
        width = 10, height = 15, dpi = 300
    )

    # split.by = "split_param"で分ける
    Feature_etiology <- FeaturePlot(
        merged_object,
        features = total_genes[i],
        pt.size = 0.5, label = TRUE, order = TRUE,
        split.by = "split_param"
    )

    Vln_etiology <- VlnPlot(
        merged_object,
        features = total_genes[i],
        pt.size = 0,
        split.by = "split_param"
    )

    plot_grid(Feature_etiology, Vln_etiology, ncol = 1, nrow = 2)
    ggsave(
        paste0("plot/", folder_name, "/total_target/Feature_Vln_etiology_", total_genes[i], ".png"),
        width = 25, height = 10, dpi = 300
    )
}