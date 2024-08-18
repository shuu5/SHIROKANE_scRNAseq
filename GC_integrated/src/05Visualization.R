library(Seurat)
library(tidyverse)

# umap_objectを読み込み
umap_object <- readRDS("data/processed/umap_object.rds")

# 各meta.dataの情報を確認
umap_object
head(umap_object)
names(umap_object@meta.data)
table(umap_object@meta.data$batch)
table(umap_object@meta.data$patient)
table(umap_object@meta.data$group)
table(umap_object@meta.data$tissue)
table(umap_object@meta.data$Diagnosis)
table(umap_object@meta.data$Lauren)
table(umap_object@meta.data$H.pylori)
table(umap_object@meta.data$celltype)
table(umap_object@meta.data$celltype, umap_object@meta.data$batch)
table(umap_object@meta.data$celltype, umap_object@meta.data$patient)
table(umap_object@meta.data$celltype, umap_object@meta.data$group)
table(umap_object@meta.data$celltype, umap_object@meta.data$tissue)
table(umap_object@meta.data$celltype, umap_object@meta.data$Diagnosis)
table(umap_object@meta.data$celltype, umap_object@meta.data$Lauren)
table(umap_object@meta.data$celltype, umap_object@meta.data$H.pylori)

### celltypeとgroupでgroup byして平均発現値を計算したseurat objectを作成
agg_group_umap_object <- AggregateExpression(umap_object, group.by = c("celltype", "group"), return.seurat = TRUE)
saveRDS(agg_group_umap_object, "data/processed/agg_group_umap_object.rds")

### celltypeとtissueでgroup byして平均発現値を計算したseurat objectを作成
agg_tissue_umap_object <- AggregateExpression(umap_object, group.by = c("celltype", "tissue"), return.seurat = TRUE)
saveRDS(agg_tissue_umap_object, "data/processed/agg_tissue_umap_object.rds")


### UMAP
# celltypeで色分け
DimPlot(
    umap_object,
    reduction = "umap_harmony",
    group.by = "celltype",
    label = TRUE,
    pt.size = 0.5
)
ggsave("output/plot/umap/umap_celltype.png", width = 10, height = 10)

# celltypeで色分け、groupでsplit
DimPlot(
    umap_object,
    reduction = "umap_harmony",
    group.by = "celltype",
    split.by = "group",
    label = TRUE,
    pt.size = 0.5
)
ggsave("output/plot/umap/umap_celltype_group.png", width = 15, height = 10)

# celltypeで色分け、tissueでsplit
DimPlot(
    umap_object,
    reduction = "umap_harmony",
    group.by = "celltype",
    split.by = "tissue",
    label = TRUE,
    pt.size = 0.5
)
ggsave("output/plot/umap/umap_celltype_tissue.png", width = 15, height = 10)

# celltypeで色分け、Diagnosisでsplit
DimPlot(
    umap_object,
    reduction = "umap_harmony",
    group.by = "celltype",
    split.by = "Diagnosis",
    label = TRUE,
    pt.size = 0.5
)
ggsave("output/plot/umap/umap_celltype_diagnosis.png", width = 25, height = 10)

# celltypeで色分け、Laurenでsplit
DimPlot(
    umap_object,
    reduction = "umap_harmony",
    group.by = "celltype",
    split.by = "Lauren",
    label = TRUE,
    pt.size = 0.5
)
ggsave("output/plot/umap/umap_celltype_lauren.png", width = 25, height = 10)

# celltypeで色分け、H.pyloriでsplit
DimPlot(
    umap_object,
    reduction = "umap_harmony",
    group.by = "celltype",
    split.by = "H.pylori",
    label = TRUE,
    pt.size = 0.5
)
ggsave("output/plot/umap/umap_celltype_hpylori.png", width = 25, height = 10)




### target geneに対するplot
tar_genes <- c("ARID5A", "ARID5B", "KRT13", "TGM2", "AXL")


### Dotplot
# celltypeで分ける
DotPlot(
    umap_object,
    features = tar_genes,
    group.by = "celltype"
) + 
theme_classic() +
theme(
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 12),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
    )
ggsave("output/plot/dotplot/dotplot_target_genes_celltype.png", width = 10, height = 10)

# celltypeごとにgroupで分ける
DotPlot(
    umap_object,
    features = tar_genes,
    group.by = "celltype",
    split.by = "group",
    cols = c("blue", "blue")
) + 
theme_classic() +
theme(
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 12),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white", colour = "white"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
    )
ggsave("output/plot/dotplot/dotplot_target_genes_celltype_group.png", width = 10, height = 10)

# celltypeごとにtissueで分ける
DotPlot(
    umap_object,
    features = tar_genes,
    group.by = "celltype",
    split.by = "tissue",
    cols = c("blue", "blue", "blue")
) +
theme_classic() +
theme(
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 12),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white", colour = "white"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
)
ggsave("output/plot/dotplot/dotplot_target_genes_celltype_tissue.png", width = 10, height = 10)


### FeaturePlot, VlnPlot
# 関数を定義
plot_gene_expression <- function(gene) {
    # 遺伝子がデータに存在するか確認
    if (gene %in% rownames(umap_object@assays$SCT$data)) {
        # FeaturePlot
        FeaturePlot(
            umap_object,
            features = gene
        )
        ggsave(paste0("output/plot/feature/featureplot_", gene, ".png"), width = 15, height = 10)

        # groupでsplit
        FeaturePlot(
            umap_object,
            features = gene,
            split.by = "group"
        )
        ggsave(paste0("output/plot/feature/featureplot_", gene, "_group.png"), width = 15, height = 10)

        # tissueでsplit
        FeaturePlot(
            umap_object,
            features = gene,
            split.by = "tissue"
        )
        ggsave(paste0("output/plot/feature/featureplot_", gene, "_tissue.png"), width = 15, height = 10)

        # Violin plot
        VlnPlot(
            umap_object,
            features = gene
        )
        ggsave(paste0("output/plot/violin/violinplot_", gene, ".png"), width = 15, height = 10)

        # groupでsplit
        VlnPlot(
            umap_object,
            features = gene,
            split.by = "group"
        )
        ggsave(paste0("output/plot/violin/violinplot_", gene, "_group.png"), width = 15, height = 10)

        # tissueでsplit
        VlnPlot(
            umap_object,
            features = gene,
            split.by = "tissue"
        )
        ggsave(paste0("output/plot/violin/violinplot_", gene, "_tissue.png"), width = 15, height = 10)
    } else {
        warning(paste("遺伝子が見つかりません:", gene))
    }
}

# tar_genesの各遺伝子についてplotを実行
for (gene in tar_genes) {
    plot_gene_expression(gene)
}


### Bar plot (平均発現値で)
# groupで分ける
for (gene in tar_genes) {
    # agg_umap_objectから遺伝子のデータを抜き出し
    if (gene %in% rownames(agg_group_umap_object@assays$SCT$data)) {
        gene_data <- as.data.frame(agg_group_umap_object@assays$SCT$data[gene, , drop = FALSE]) %>%
            rownames_to_column(var = "gene") %>%
            pivot_longer(cols = -gene, names_to = "celltype_group", values_to = "value")

        # celltypeとgroupを分割
        long_data <- gene_data %>%
            separate(celltype_group, into = c("celltype", "group"), sep = "_")

        # Bar plotの作成
        bar_plot <- ggplot(long_data, aes(x = celltype, y = value, fill = group)) +
            geom_bar(stat = "identity", position = "dodge", color = "black") +  # 外枠を黒線に設定
            scale_fill_manual(values = c("Non-Tumor" = "gray", "Tumor" = "red")) +
            theme_classic() +
            labs(title = paste0(gene, " Average Expression"), x = "Cell Type", y = "Average Expression") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        ggsave(paste0("output/plot/barplot/barplot_", gene, "_celltype_group.png"), plot = bar_plot, width = 10, height = 10)
    } else {
        warning(paste("遺伝子が見つかりません:", gene))
    }
}

# tissueで分ける
for (gene in tar_genes) {
    # agg_tissue_umap_objectから遺伝子のデータを抜き出し
    if (gene %in% rownames(agg_tissue_umap_object@assays$SCT$data)) {
        gene_data <- as.data.frame(agg_tissue_umap_object@assays$SCT$data[gene, , drop = FALSE]) %>%
            rownames_to_column(var = "gene") %>%
            pivot_longer(cols = -gene, names_to = "celltype_tissue", values_to = "value")
    
        # celltypeとtissueを分割
        long_data <- gene_data %>%
            separate(celltype_tissue, into = c("celltype", "tissue"), sep = "_")
    
        # Bar plotの作成
        bar_plot <- ggplot(long_data, aes(x = celltype, y = value, fill = tissue)) +
            geom_bar(stat = "identity", position = "dodge", color = "black") +  # 外枠を黒線に設定
            scale_fill_manual(values = c("Normal" = "gray", "Paratumor" = "darkgreen", "Tumor" = "red")) +
            theme_classic() +
            labs(title = paste0(gene, " Average Expression"), x = "Cell Type", y = "Average Expression") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
        ggsave(paste0("output/plot/barplot/barplot_", gene, "_celltype_tissue.png"), plot = bar_plot, width = 10, height = 10)
    } else {
        warning(paste("遺伝子が見つかりません:", gene))
    }
}