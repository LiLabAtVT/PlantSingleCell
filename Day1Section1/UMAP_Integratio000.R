library(Seurat)
library(ggplot2)

# Set data directory and sample names
data_dir <- "/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Razan_Nina_Integration_Analysis/Data"
sample_files <- c("Neg_hpi1.rds", "Neg_hpi2.rds", "Pos_hpi1.rds", "Pos_hpi2.rds")
resolutions <- seq(0.2, 0.9, 0.1)

# Function to process a Seurat object
process_seurat_object <- function(obj_file, data_dir, resolutions) {
  obj_path <- file.path(data_dir, obj_file)
  seurat_obj <- readRDS(obj_path)
  
  # Run PCA if not already present
  if (is.null(seurat_obj@reductions$pca)) {
    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  }
  
  # Run UMAP if not already present
  if (is.null(seurat_obj@reductions$umap)) {
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
  }
  
  # Prepare PDF for UMAP plots
  pdf(file = file.path(data_dir, paste0(tools::file_path_sans_ext(obj_file), "_UMAPs.pdf")), width = 8, height = 6)
  for (res in resolutions) {
    seurat_obj <- FindClusters(seurat_obj, resolution = res)
    p <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, label.size = 5) +
      ggtitle(paste0(obj_file, " | Resolution: ", res)) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
  }
  dev.off()
  
  # Save processed object
  saveRDS(seurat_obj, file = file.path(data_dir, paste0(tools::file_path_sans_ext(obj_file), "_clustered0.2_0.9.rds")))
}

# Process all objects
for (obj in obj_names) {
  process_seurat_object(obj, data_dir, resolutions)
}
