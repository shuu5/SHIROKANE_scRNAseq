library(Seurat)
library(ComplexHeatmap)
library(clustermole)
library(tidyverse)

# data/processed/subclusters/内のsubclusteringしたオブジェクトを読み込み
subclust_list <- lapply(list.files("data/processed/subclusters/", pattern = "*.rds", full.names = TRUE), readRDS)
names(subclust_list) <- gsub("data/processed/subclusters/(.*)_subclust.rds", "\\1", list.files("data/processed/subclusters/", pattern = "*.rds"))
names(subclust_list) <- gsub("_subclust.rds", "", names(subclust_list))

# subclust_listの各オブジェクトに対して全クラスタの変動遺伝子を抽出
cluster_markers_list <- list()
for (name in names(subclust_list)) {
  umap_object <- subclust_list[[name]]
  print(paste0("FindAllMarkers: ", name))
  umap_object <- PrepSCTFindMarkers(umap_object)
  cluster_markers <- FindAllMarkers(
      umap_object,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.5
  )
  cluster_markers_list[[name]] <- cluster_markers

  # markerをrdsとcsvに保存
  saveRDS(cluster_markers, file = paste0("output/rds/annotation/subclusters/", name, "_cluster_markers.rds"))
  write.csv(cluster_markers, file = paste0("output/csv/annotation/subclusters/", name, "_cluster_markers.csv"))

  # 各クラスターのDEGの数を確認
  print(paste0("cluster_DEG_number"))
  print(table(cluster_markers$cluster))

  # クラスターごとにDEGを分割してリストにまとめる
  markers_list <- split(cluster_markers$gene, cluster_markers$cluster)

  # リスト内要素の共通項のパターンを集計
  common_pattern <- ComplexHeatmap::make_comb_mat(markers_list)

  ### UpSet Plot
  print(paste0("processing UpSet"))
  png(paste0("output/plot/UpSet/subclusters/", name, "_UpSet_plot.png"))
  ComplexHeatmap::UpSet(common_pattern)
  dev.off()

  # size順に並び変え
  png(paste0("output/plot/UpSet/subclusters/", name, "_UpSet_plot_size.png"))
  ComplexHeatmap::UpSet(
      common_pattern, 
      comb_order = order(ComplexHeatmap::comb_size(common_pattern))
  )
  dev.off()

  # クラスター名順に並び変え
  png(paste0("output/plot/UpSet/subclusters/", name, "_UpSet_plot_cluster.png"))
  ComplexHeatmap::UpSet(
      common_pattern,
      set_order = sort(names(markers_list))
  )   
  dev.off()

  ### 各遺伝子とパターンコード、他との重複数をデータフレームにまとめる
  print(paste0("processing extract_comb: ", name))
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
  write.csv(genes_df, file = paste0("output/csv/annotation/subclusters/", name, "_genes_df.csv"))
  saveRDS(genes_df, file = paste0("output/rds/annotation/subclusters/", name, "_genes_df.rds"))
}


