library(Seurat)
library(ComplexHeatmap)
library(clustermole)
library(tidyverse)

# umap_objectを読み込み
umap_object <- readRDS("data/processed/umap_object.rds")

# 全クラスタの変動遺伝子を抽出
umap_object <- PrepSCTFindMarkers(umap_object)
cluster_markers <- FindAllMarkers(
    umap_object,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.5
)

# markerをrdsとcsvに保存
saveRDS(cluster_markers, file = "output/rds/annotation/cluster_markers.rds")
write.csv(cluster_markers, file = "output/csv/annotation/cluster_markers.csv")

# 各クラスターのDEGの数を確認
table(cluster_markers$cluster)

# クラスターごとにDEGを分割してリストにまとめる
markers_list <- split(cluster_markers$gene, cluster_markers$cluster)

# リスト内要素の共通項のパターンを集計
common_pattern <- ComplexHeatmap::make_comb_mat(markers_list)



### UpSet Plot
png("output/plot/UpSet/UpSet_plot.png")
ComplexHeatmap::UpSet(common_pattern)
dev.off()

# size順に並び変え
png("output/plot/UpSet/UpSet_plot_size.png")
ComplexHeatmap::UpSet(
    common_pattern, 
    comb_order = order(ComplexHeatmap::comb_size(common_pattern))
    )
dev.off()

# クラスター名順に並び変え
png("output/plot/UpSet/UpSet_plot_cluster.png")
ComplexHeatmap::UpSet(
    common_pattern,
    set_order = sort(names(markers_list))
    )   
dev.off()


### 各遺伝子とパターンコード、他との重複数をデータフレームにまとめる
# パターンコード一覧を取得
patterns <- comb_name(m = common_pattern)

# 空のリストを用意。名前も付けておく
extract_genes <- vector(mode = "list", length = length(patterns))
names(extract_genes) <- patterns

# パターンコードごとにextract_combで取り出して、リストに保存
for(i in patterns){
  extract_genes[[i]] <- extract_comb(m = common_pattern, comb_name = i)
}

# sapplyでリストにlength関数を繰り返し処理した結果を、orderで降順に並び替え
extract_genes <- extract_genes[order(sapply(extract_genes, length), decreasing = T)]

# 遺伝子名とパターンコードについてのデータフレームを作成する機能をリストに繰り返し処理
tmp <- mapply(function(x, y){
  df <- data.frame(genes = x,  # 遺伝子ベクトルが入ったリストから、遺伝子ベクトルを取り出し
                   code = rep(y, times = length(x)) # リストの名前 (パターンコード)を遺伝子ベクトルの要素数だけ繰り返し
  )
  return(df)
  },
  x = extract_genes,
  y = names(extract_genes),
  SIMPLIFY = F # 結果をリストとして出力させる
  )

# リストの要素をrbindで縦につなげる
genes_df <- do.call(rbind, tmp)
rownames(genes_df) <- genes_df$genes

# str_counts関数でパターンコード中の1の数を調べる。
genes_df$num_overlap <- stringr::str_count(string = genes_df$code, pattern = "1")

# genes_dfをcsvとrdsに保存
write.csv(genes_df, file = "output/csv/annotation/genes_df.csv")
saveRDS(genes_df, file = "output/rds/annotation/genes_df.rds")


### clustemoleでannotattion
# 各クラスターについてoverlapしているセルタイプを取得
res_list <- lapply(unique(cluster_markers$cluster), function(x) {
    genes <- cluster_markers[cluster_markers$cluster == x, "gene"]
    res <- clustermole_overlaps(genes = genes, species = "hs")
    return(res)
})
names(res_list) <- unique(cluster_markers$cluster)

# 結果の確認と上位10個のcelltypeを抽出
top_celltypes <- lapply(res_list, function(x) {
    head(x$celltype, 10)
})

# 上位10個のcelltypeをデータフレームに変換
top_celltypes_df <- do.call(rbind, lapply(names(top_celltypes), function(cluster) {
    data.frame(cluster = cluster, celltype = top_celltypes[[cluster]])
}))

# CSV形式で保存
write.csv(top_celltypes_df, file = "output/csv/annotation/top_celltypes_per_cluster.csv", row.names = FALSE)



### clusterの平均発現値でannotationする方法
aveexp <- AverageExpression(umap_object, layer = "counts")$SCT
aveexp <- log1p(aveexp)
aveexp <- as.matrix(aveexp)

enrich_res <- clustermole_enrichment(expr_mat = aveexp, species = "hs")

# 平均発現値の上位10個のcelltypeを抽出
top_celltypes_aveexp <- lapply(unique(enrich_res$cluster), function(x) {
    # 修正: 上位10個のcelltypeをリストとして取得
    top_celltypes <- head(enrich_res[enrich_res$cluster == x, "celltype"], 10)
    # 修正: clusterとcelltypeをデータフレームに変換
    data.frame(cluster = x, celltype = top_celltypes)
})

# 上位10個のcelltypeをデータフレームに変換
top_celltypes_aveexp_df <- do.call(rbind, top_celltypes_aveexp)

# CSV形式で保存
write.csv(top_celltypes_aveexp_df, file = "output/csv/annotation/top_celltypes_aveexp_per_cluster.csv", row.names = FALSE)