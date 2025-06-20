library(Seurat)
library(pheatmap)
library(qs)

set.seed(123456)
options(future.globals.maxSize = 200 * 1024^3)  # Set to 200GB

## Function definitions
heatmap_cluster <- function(corr_mat, k, mat_type, layer) {
  # save heatmap
  dir_name <- file.path("./Root/Correlation heatmaps", mat_type)
  print(dir_name)
  print(file.path(dir_name, paste0("marker_corr_", layer, "_", k, ".png")))

  if (!dir.exists(dir_name)) {
        dir.create(dir_name, recursive = TRUE)
  }
  output_path <- file.path(dir_name, paste0("marker_corr_", layer, "_", k, ".png"))
  png(output_path, width = 600, height = 600, type = "cairo")
  corr_pheatmap_res <- pheatmap(corr_mat, cutree_rows = k, cutree_cols = k)
  dev.off()

   # Extract cluster assignments
  row_clusters <- cutree(corr_pheatmap_res$tree_row, k = k)

  # Save clusters to CSV
  dir_name <- file.path("./Root/Tree cut clusters CSVs", mat_type)
  if (!dir.exists(dir_name)) {
        dir.create(dir_name, recursive = TRUE)
  }
  output_path <- file.path(dir_name, paste0("marker_clusters_", layer, "_", k, ".csv"))
  write.csv(data.frame(Gene = names(row_clusters), Cluster = row_clusters),
   output_path, row.names = FALSE)

  return(corr_pheatmap_res)
}

save_corr_result_ordered <- function(corr_mat, row_order, col_order, mat_type, layer){
  corr_mat_ordered <- corr_mat[row_order, col_order]

  dir_name <- file.path("./Root/Correlation CSVs/", mat_type)
  if (!dir.exists(dir_name)) {
        dir.create(dir_name, recursive = TRUE)
  }

  output_path <- file.path(dir_name, paste0("/marker_corr_", layer,"_ordered.csv"))
  write.csv(corr_mat_ordered, output_path)
}



#### marker genes
print("Marker genes.......")
genes_of_interest <- read.csv("/projects/songli_lab/PlantSingleCell2025/Day_1/Session_2/Inputs/Root genes.csv")
gene_ids_of_interest <- genes_of_interest$GeneID
print("# of genes of interest")
print(length(gene_ids_of_interest))

# read seurat objects
integrated_cca <- qread("./Root/integrated_cca.qs")
integrated_rpca <- qread("./Root/integrated_rpca.qs")
merged <- qread("./Root/merged.qs")

# common genes between all genes and genes of interest
genes <- as.list(rownames(merged))
gene_ids_of_interest <-  as.list(gene_ids_of_interest)
common_genes_with_all <- unlist(intersect(genes, gene_ids_of_interest))

print("# common genes between merged and genes of interest:")
print(length(common_genes_with_all))
missing_genes <- setdiff(gene_ids_of_interest, common_genes_with_all)
print("# missing genes:")
print(length(missing_genes))

# common genes between integrated genes and genes of interest
genes <- as.list(rownames(integrated_cca))
gene_ids_of_interest <-  as.list(gene_ids_of_interest)
common_genes_with_integrated <- unlist(intersect(genes, gene_ids_of_interest))
print("# common genes between integrated genes of interest:")
print(length(common_genes_with_integrated))
missing_genes <- setdiff(gene_ids_of_interest, common_genes_with_integrated)
print("# of missing genes:")
print(length(missing_genes))


#### Carrelation for common genes
### Merged
merged <- NormalizeData(merged, verbose = TRUE)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
merged <- ScaleData(merged, features = rownames(merged), verbose = TRUE)

row_indices <- which(rownames(merged) %in% common_genes_with_integrated)
marker_merged <- merged[row_indices, ]

## Merged: Count layer
marker_merged_mat_count <- as.matrix(GetAssayData(marker_merged, layer = "count"))
marker_merged_corr_count <- cor(t(log2(marker_merged_mat_count+1)))
# save heatmap plots
corr_pheatmap_res <- heatmap_cluster(marker_merged_corr_count, 4, "merged", "count")
# get correlation data based on clustered order
row_order <- rownames(marker_merged_corr_count)[corr_pheatmap_res$tree_row$order]
col_order <- colnames(marker_merged_corr_count)[corr_pheatmap_res$tree_col$order]

heatmap_cluster(marker_merged_corr_count, 6, "merged", "count")
heatmap_cluster(marker_merged_corr_count, 8, "merged", "count")

# save ordered correlation reults
save_corr_result_ordered(marker_merged_corr_count, row_order, col_order, "merged", "count")
# -----------

## Merged: Data layer
marker_merged_mat_data <- as.matrix(GetAssayData(marker_merged, layer = "data"))
marker_merged_corr_data <- cor(t(marker_merged_mat_data))
print("dim of correlation matrix for Merged (data):")
print(dim(marker_merged_corr_data))
# save heatmap plots
heatmap_cluster(marker_merged_corr_data, 4, "merged", "data")
heatmap_cluster(marker_merged_corr_data, 6, "merged", "data")
heatmap_cluster(marker_merged_corr_data, 8, "merged", "data")
# save ordered correlation reults
save_corr_result_ordered(marker_merged_corr_data, row_order, col_order, "merged", "data")
# -----------

## Merged: Scale.data layer
marker_merged_mat_scale_data <- as.matrix(GetAssayData(marker_merged, layer = "scale.data"))
marker_merged_corr_scale_data <- cor(t(marker_merged_mat_scale_data))
print("dim of correlation matrix for Merged (scale.data):")
print(dim(marker_merged_corr_scale_data))
# save heatmap plots
heatmap_cluster(marker_merged_corr_scale_data, 4, "merged", "scale.data")
heatmap_cluster(marker_merged_corr_scale_data, 6, "merged", "scale.data")
heatmap_cluster(marker_merged_corr_scale_data, 8, "merged", "scale.data")
# save ordered correlation reults
save_corr_result_ordered(marker_merged_corr_scale_data, row_order, col_order, "merged", "scale.data")
# -----------

### Integrated - CCA
row_indices_cca <- which(rownames(integrated_cca) %in% common_genes_with_integrated)
marker_integrated_cca <- integrated_cca[row_indices_cca, ]

## Integrated - CCA: Scale.data layer
marker_integrated_cca_mat_scale_data <- as.matrix(GetAssayData(marker_integrated_cca, assay = "integrated", layer = "scale.data"))
marker_integrated_cca_corr_scale_data <- cor(t(marker_integrated_cca_mat_scale_data))
print("dim of correlation matrix for Integrated_CCA (scale.data):")
print(dim(marker_integrated_cca_corr_scale_data))
# save heatmap plots
heatmap_cluster(marker_integrated_cca_corr_scale_data, 4, "integrated_CCA", "scale.data")
heatmap_cluster(marker_integrated_cca_corr_scale_data, 6, "integrated_CCA", "scale.data")
heatmap_cluster(marker_integrated_cca_corr_scale_data, 8, "integrated_CCA", "scale.data")
# save ordered correlation reults
save_corr_result_ordered(marker_integrated_cca_corr_scale_data, row_order, col_order, "integrated_CCA", "scale.data")
# -----------

### Integrated - RPCA
row_indices_rpca <- which(rownames(integrated_rpca) %in% common_genes_with_integrated)
marker_integrated_rpca <- integrated_rpca[row_indices_rpca, ]

## Integrated - RPCA: Scale.data layer
marker_integrated_rpca_mat_scale_data <- as.matrix(GetAssayData(marker_integrated_rpca, assay = "integrated", layer = "scale.data"))
marker_integrated_rpca_corr_scale_data <- cor(t(marker_integrated_rpca_mat_scale_data))
print("dim of correlation matrix for Integrated_RPCA (scale.data):")
print(dim(marker_integrated_rpca_corr_scale_data))
# save heatmap plots
heatmap_cluster(marker_integrated_rpca_corr_scale_data, 4, "integrated_RPCA", "scale.data")
heatmap_cluster(marker_integrated_rpca_corr_scale_data, 6, "integrated_RPCA", "scale.data")
heatmap_cluster(marker_integrated_rpca_corr_scale_data, 8, "integrated_RPCA", "scale.data")
# save ordered correlation reults
save_corr_result_ordered(marker_integrated_rpca_corr_scale_data, row_order, col_order, "integrated_RPCA", "scale.data")