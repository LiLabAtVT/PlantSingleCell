# Integration of Pathogen and Non-pathogen scRNA-seq data

library(Seurat)
library(sctransform)

# Load raw 10X count matrices
Neg1 <- Read10X(data.dir = "/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/CellRanger_Output_Mapped_Multi_Ref/PLant_Pathogen_Mapping_Neg24hpi_1_scRNAseq_Multi_Ref/outs/filtered_feature_bc_matrix")
Pos1 <- Read10X(data.dir = "/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/Convert_BAM_to_FASTQ/Mapping_Pos24hpi_1_scRNAseq_ATH_Read_Only/outs/filtered_feature_bc_matrix")

# Clean gene names
rownames(Neg1) <- gsub("Arabidopsis_thaliana_gene:", "", rownames(Neg1))
rownames(Pos1) <- gsub("gene:", "", rownames(Pos1))

# Seurat object processing function
process_seurat_object <- function(counts_data, project_name,
                                  feature_upper, nfeatures,
                                  pca_dims, clustering_resolution) {
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts_data, project = project_name, 
                                   min.cells = 3, min.features = 200)
  
  # Compute mitochondrial content
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^ATMG")
  
  # Violin plot before filtering
  pdf(paste0("/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/RScript/", project_name, "_QC_violin_pre_filter.pdf"))
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  dev.off()
  
  # Filter cells
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 100 & 
                         nFeature_RNA < feature_upper & 
                         percent.mt < 5)
  
  # Violin plot after filtering
  pdf(paste0("/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/RScript/", project_name, "_QC_violin_post_filter.pdf"))
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  dev.off()
  
  # Normalization and variable features
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = nfeatures)
  
  # Scale, PCA, clustering
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  seurat_obj <- FindNeighbors(seurat_obj, dims = pca_dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = clustering_resolution)
  seurat_obj <- RunUMAP(seurat_obj, dims = pca_dims)
  
  # Save object
  saveRDS(seurat_obj, paste0("/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/RScript/", project_name, ".rds"))
  
  return(seurat_obj)
}

# Process both samples
Neg24hpi_1 <- process_seurat_object(Neg1, project_name = "Neg_hpi1", 
                                    feature_upper = 6000, nfeatures = 2000,
                                    pca_dims = 1:20, clustering_resolution = 0.22)

Pos24hpi_1 <- process_seurat_object(Pos1, project_name = "Pos_hpi1", 
                                    feature_upper = 3500, nfeatures = 2000,
                                    pca_dims = 1:20, clustering_resolution = 0.22)

# Merge and integrate
crossSamples <- SplitObject(merge(Neg24hpi_1, y = c(Pos24hpi_1)), split.by = "orig.ident")

for (i in names(crossSamples)) {
  crossSamples[[i]] <- SCTransform(crossSamples[[i]], verbose = FALSE)
}

features <- SelectIntegrationFeatures(object.list = crossSamples, nfeatures = 10000)
crossSamples <- PrepSCTIntegration(object.list = crossSamples, anchor.features = features)
crossSamples.anchors <- FindIntegrationAnchors(object.list = crossSamples, normalization.method = "SCT", anchor.features = features)
crossSamples.combine <- IntegrateData(anchorset = crossSamples.anchors, normalization.method = "SCT")

# Run dimensional reduction and clustering
crossSamples.combine <- RunPCA(crossSamples.combine, verbose = FALSE)
crossSamples.combine <- RunUMAP(crossSamples.combine, reduction = "pca", dims = 1:30)
crossSamples.combine <- FindNeighbors(crossSamples.combine, reduction = "pca", dims = 1:30)
crossSamples.combine <- FindClusters(crossSamples.combine, resolution = 0.3)

# Save integrated object
saveRDS(crossSamples.combine, "/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/RScript/Patho_Nonpatho_integration.Rds")

# UMAP visualizations
pdf("/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/RScript/Patho_Nonpatho.pdf", width = 8, height = 5)
DimPlot(crossSamples.combine, reduction = "umap", group.by = "orig.ident", cols = c("lightblue", "blue", "#FFC073", "#FF8C00"))
dev.off()

pdf("/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/RScript/Patho_Nonpatho_sepColor.pdf", width = 8, height = 5)
DimPlot(crossSamples.combine, reduction = "umap", split.by = "orig.ident", group.by = "orig.ident", cols = c("lightblue", "blue", "#FFC073", "#FF8C00"))
dev.off()

pdf("/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/RScript/Patho_Nonpatho_label.pdf", width = 8, height = 5)
DimPlot(crossSamples.combine, reduction = "umap", label = TRUE)
dev.off()

pdf("/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/RScript/Patho_Nonpatho_sep.pdf", width = 8, height = 5)
DimPlot(crossSamples.combine, reduction = "umap", split.by = "orig.ident")
dev.off()
