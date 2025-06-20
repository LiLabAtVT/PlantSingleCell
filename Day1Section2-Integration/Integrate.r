library(Seurat)
library(sctransform)
library(glmGamPoi) # to speed up the SCTransform process
library(patchwork) # using wrap_plots function
library(ggplot2)
library(RColorBrewer)

set.seed(123456)
options(future.globals.maxSize = 200 * 1024^3)  # Set to 200GB

# Mapping of detailed cell‐types → biological clusters
cluster2celltype <- list(
  "Epidermal" = c("trichoblast", "atrichoblast/epidermis",  "atrichoblast",  "root cap"),
  "Ground"    = c("cortex", "endodermis (Casparian strip)", "endodermis", "suberized endodermis"),
  "Vascular"  = c("pericycle", "companion cell", "phloem parenchyma", "xylem/procambium"),
  "Meristem"  = c("dividing cell", "initial cell")
)
# invert to look up cluster by annotation
celltype2cluster <- setNames(
  rep(names(cluster2celltype), lengths(cluster2celltype)),
  unlist(cluster2celltype)
)

## Function definitions
get_color_pallette <- function(clustertype) {
  if (clustertype == "per_celltype") {
    # Set color for each cell type but make colors of each cluster close
    library(RColorBrewer)
    
    palette_names <- list(
      Epidermal = "Blues",
      Ground    = "Greens",
      Vascular  = "YlOrBr",
      Meristem  = "Greys"
    )

    group_palettes <- lapply(names(cluster2celltype), function(g){
      n <- length(cluster2celltype[[g]])
      pal <- palette_names[[g]]

      if (n >= 3) {
        full_pal <- brewer.pal(9, pal)
        indices <- seq(9, 1, by = -2)[1:n]
        full_pal[indices]
      } else {
        brewer.pal(3, pal)[(3 - n + 1):3]
      }
    })

    names(group_palettes) <- names(cluster2celltype)

    cols <- unlist(lapply(names(group_palettes), function(g){
      setNames(group_palettes[[g]], cluster2celltype[[g]])
    }), use.names = TRUE)

  } else if (clustertype == "per_cluster") {
    clusters <- names(cluster2celltype)

    cluster_cols <- c(
      brewer.pal(9, "Blues")[5],   # Epidermal → medium blue
      brewer.pal(9, "Greens")[5],  # Ground → medium green
      brewer.pal(9, "YlOrBr")[5],  # Vascular → medium yellow-orange-brown
      brewer.pal(9, "Greys")[5]    # Meristem → medium grey
    )


    names(cluster_cols) <- clusters

    cols <- cluster_cols[celltype2cluster]
    celltypes <- names(celltype2cluster)
    names(cols) <- celltypes

  }
  else {
    stop("Unknown clustertype. Choose from 'per_celltype', 'per_cluster', or 'per_custom_cluster'.")
  }

  return(cols)
}

dimplot <- function(obj, 
                               group_by,
                               split_by = NULL,
                               pallette = NULL, 
                               sample_name = NULL, 
                               title_prefix = NULL,
                               file_path = NULL,
                               width = 14, height = 6) {
  
  # title
  plot_title <- if (is.null(title_prefix)) {
    sample_name
  } else {
    paste(title_prefix, sample_name)
  }
  
  plot_args <- list(
    object     = obj,
    reduction  = "umap",
    group.by   = group_by,
    pt.size    = 0.5,
    label      = FALSE,
    label.size = 3
  )
  
  # add split.by if provided
  if (!is.null(split_by)) {
    plot_args$split.by <- split_by
  }
  # add cols if provided
  if (!is.null(pallette)) {
    plot_args$cols <- pallette
  }

  p <- do.call(DimPlot, plot_args) +
    theme_minimal() +
    ggtitle(plot_title)

  # save if file_path is not null
  if (!is.null(file_path)) {
    ggsave(
        filename = file_path,
        plot = p,
        width = width,
        height = height,
        dpi = 300,
        units = "in"
        )
      return(invisible(TRUE))
  }
    
  return(p)
}

dimplot_combined <- function(
  obj_list,
  group_by         = "integrated_annotation",
  title_prefix     = NULL,
  pallette         = NULL,
  file             = NULL,
  nrow  = 1, ncol = length(obj_list),
  width  = 14, height = 6, dpi = 300
) {

  if (is.null(names(obj_list))) {
    stop("`obj_list` must be a named list (names will be used as plot titles).")
  }
  
  plot_list <- lapply(seq_along(obj_list), function(i) {
    sample_name <- names(obj_list)[i]
    obj <- obj_list[[i]]
    obj <- SetIdent(obj, value = group_by)

    p <- dimplot(
      obj          = obj,
      group_by     = group_by,
      pallette     = pallette,
      sample_name  = sample_name,
      title_prefix = title_prefix
      )

      if (i > 1) {
        p <- p + theme(legend.position = "none") + guides(color = "none")
      }

      p

    })
    
    combined <- wrap_plots(plot_list, ncol = ncol, nrow = nrow)
    
    if (!is.null(file)) {

      if (!dir.exists(dirname(file))) {
        dir.create(dirname(file), recursive = TRUE)
      }

      ggsave(
        filename = file,
        plot = combined,
        width = width,
        height = height,
        dpi = dpi,
        units = "in"
        )
      return(invisible(TRUE))
    }
    return(combined)
  }


### 1. Loading data
rds_files <- c("/projects/songli_lab/PlantSingleCell2025/Day_1/Session_2/Inputs/DAG6_root.slim.rds", 
               "/projects/songli_lab/PlantSingleCell2025/Day_1/Session_2/Inputs/DAG11_root.slim.rds")

seurat_list <- lapply(rds_files, readRDS)

print("Original seurat object dimensions:")
print(sapply(seurat_list, dim))

print("Original seurat object assays:")
print(sapply(seurat_list, function(x) names(x@assays)))

# assign unique orig.ident based on file names or indices
orig_ident_map <- c("Root (6 DAG)", "Root (11 DAG)")

for (i in seq_along(seurat_list)) {
  seurat_list[[i]]$orig.ident <- orig_ident_map[i]
}
print("unique(seurat_object$orig.ident)")
print(sapply(seurat_list, function(x) unique(x$orig.ident)))

### 2. Merge-then-split: to have a consistent gene set across all datasets.
  #  merge the first Seurat object with the rest of the Seurat objects
  #  Cells: All cells from all datasets are retained.
  #  Gene Expression: If a gene was present in a particular dataset, its expression value for that cell will remain unchanged.
  #  Missing Genes: If a gene was not measured in a particular dataset, the expression value for that gene in the cells from that dataset will be set to NA (.).

#  each Seurat object will have its 'origin.ident' added to its cell IDs during the merge, preventing conflicts between cell names
merged <- merge(seurat_list[[1]], y = seurat_list[-1], 
               add.cell.ids = sapply(seurat_list, function(x) unique(x$orig.ident)))

if (!dir.exists("./Root")) {
  dir.create("./Root")
}
saveRDS(merged, file = "./Root/merged.rds")

print("Merged seurat object dimensions:")
print(dim(merged)) # 22403 genes are common between these two samples

## split merged into multiple smaller Seurat objects based on 'orig.ident'
# gives same cell numebrs as original samples but integrates all genes
cross_samples <- SplitObject(merged, split.by = "orig.ident")
print("Sample dimensions after merge-split:")
print(sapply(cross_samples, dim))

cross_samples <- lapply(cross_samples, function(obj) {
  # set cell type as ident for each sample
  obj <- SetIdent(obj, value = "integrated_annotation")
  
  # data normalization and dimension reduction
  obj <- NormalizeData(obj, verbose = TRUE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, features = rownames(obj), verbose = TRUE)
  obj <- RunPCA(obj, verbose = TRUE)
  obj <- RunUMAP(obj, dims = 1:30, verbose = TRUE)

  # add the new meta column "celltype_cluster"
  fine <- as.character(Idents(obj))
  cells <- Cells(obj)
  names(fine) <- cells
  cluster <- celltype2cluster[fine]
  names(cluster) <- cells
  obj$celltype_cluster <- cluster

  print(unique(obj@meta.data$celltype_cluster))
  
  return(obj)
  })

# 3. SCTransform normalization and variance stabilization on each sample and set default assay to "SCT"
  # variance stabilization: of the gene expression values to minimize the impact of highly variable genes, making downstream analyses like clustering and PCA more reliable. 
  # normalization (using a method based on regularized negative binomial regression for each gene) to account for technical variations and differences in sequencing depth across cells.)
  # It also applies two filters that remove genes from the output assay:
  # 1. Low‐expression filtering via min_cells: discards any gene that is detected in fewer than 5 cells (min_cells = 5)
  # 2. Only returns the set of variable features (the genes whose residuals show high variance across cells) in the assay slots (return.only.var.genes = TRUE)
sc_cross_samples <- lapply(cross_samples, function(obj) {
  obj <- SCTransform(obj, verbose=TRUE)

  # PCA will on the  variance-stabilized data from the SCT assay
  obj <- RunPCA(obj, verbose = TRUE)
  # UMAP on the PCA results from the SCT assay
  obj <- RunUMAP(obj, dims = 1:30, verbose = TRUE)
  return(obj)
})

print("Sample dimensions after merge-split and SCT:")
print(sapply(sc_cross_samples, dim))
print("Sample assays after merge-split and SCT:")
print(sapply(sc_cross_samples, function(x) names(x@assays)))
print("Default assay after merge-split and SCT:")
sapply(sc_cross_samples, DefaultAssay)

print("Plot UMAP before integration (on the SCT-normalized objects)")
dimplot_combined(
  obj_list    = sc_cross_samples,
  group_by    = "integrated_annotation",
  pallette    = get_color_pallette("per_celltype"),
  file        = "./Root/UMAP plots/umap_individual_samples_per_celltype.png",
)
dimplot_combined(
  obj_list    = sc_cross_samples,
  group_by    = "integrated_annotation",
  pallette    = get_color_pallette("per_cluster"),
  file        = "./Root/UMAP plots/umap_individual_samples_per_cluster.png",
)

# 4. Select integration features which identifies shared, HVGs across samples
  # to help minimize the influence of batch effects (unwanted technical differences between datasets) that can obscure biological signals
  # 4.1. selects the top 'nfeatures' HVG within each individual sample
  # 4.2. finds genes that are variable in most samples, not necessarili all (intersection of HVGs)
features <- SelectIntegrationFeatures(sc_cross_samples, nfeatures = 5000)
print("length(IntegrationFeatures) with nfeatures>5000 - SCT:")
print(length(features))

# 5. Recalculate SCT residuals only for the genes provided in 'anchor.features' and store residuals in the SCT assay slot inside each object
# output: a prepared list of Seurat objects ready for anchor finding and integration.
sc_cross_samples <- PrepSCTIntegration(sc_cross_samples, anchor.features = features)

# 6. Find integration anchors
# from papeR: identifies "anchors", which are correspondences between cells in different experiments
# that represent shared biological states. These anchors facilitate harmonizing datasets into a single
# reference and transferring information like cell labels and gene expression patterns between different types of single-cell measurements. 

# finding pairs of biologically similar cells across samples (cells that match each other across samples) based on selected features (genes)
  # with anchors, Seurat can preserve biological signals and correct for technical noise.
# reduction = CCA (Canonical Correlation Analysis)
sc_cross_samples.cca_anchors <- FindIntegrationAnchors(
  object.list = sc_cross_samples,
  normalization.method = "SCT", # use SCT residuals for integration (not log-normalized counts)
  anchor.features = features, # limits the search to the selected integration features
  reduction = "cca",
  verbose = TRUE
)
print("sc_cross_samples.CCA_anchors:")
print(sc_cross_samples.cca_anchors)
print(head(sc_cross_samples.cca_anchors@anchors))

# reduction = RPCA (Reciprocal PCA) - fast
sc_cross_samples.rpca_anchors <- FindIntegrationAnchors(
  object.list = sc_cross_samples,
  normalization.method = "SCT", 
  anchor.features = features,
  reduction = "rpca",
  verbose = TRUE
)
print("sc_cross_samples.RPCA_anchors:")
print(sc_cross_samples.rpca_anchors)
print(head(sc_cross_samples.rpca_anchors@anchors))

# 7. Integrate data
# using anchors to mathematically correct batch effects and align datasets into a common space.
integrated_cca <- IntegrateData(
  anchorset = sc_cross_samples.cca_anchors,
  normalization.method = "SCT"
)
print("Assays(integrated_cca):")
print(Assays(integrated_cca))

print("DefaultAssay(integrated_cca):")
print(DefaultAssay(integrated_cca))

print("Dimension of integrated assay (CCA)")
print(dim(integrated_cca@assays$integrated))

integrated_rpca <- IntegrateData(
  anchorset = sc_cross_samples.rpca_anchors,
  normalization.method = "SCT"
)

# PCA and UMAP on the default assay ("integrated")
integrated_cca <- RunPCA(integrated_cca, verbose = TRUE)
integrated_cca <- RunUMAP(integrated_cca, reduction = "pca", dims = 1:30)
saveRDS(integrated_cca, "./Root/integrated_cca.rds")

integrated_rpca <- RunPCA(integrated_rpca, verbose = TRUE)
integrated_rpca <- RunUMAP(integrated_rpca, reduction = "pca", dims = 1:30)
saveRDS(integrated_rpca, "./Root/integrated_rpca.rds")

# UMAP plots
dimplot(
  obj       = integrated_cca,
  group_by  = "orig.ident",
  file_path = "./Root/UMAP plots/integrated_CCA_UMAP_per-sample.png",
  width = 7, height = 5
  )
dimplot(
  obj       = integrated_cca,
  group_by  = "integrated_annotation",
  split_by  = "orig.ident",
  pallette  = get_color_pallette("per_celltype"),
  file_path = "./Root/UMAP plots/split_integrated_CCA_UMAP_per_celltype.png"
  )
dimplot(
  obj       = integrated_cca,
  group_by  = "integrated_annotation",
  split_by  = "orig.ident",
  pallette  = get_color_pallette("per_cluster"),
  file_path = "./Root/UMAP plots/split_integrated_CCA_UMAP_per_cluster.png"
  )


  dimplot(
  obj       = integrated_rpca,
  group_by  = "orig.ident",
  file_path = "./Root/UMAP plots/integrated_RPCA_UMAP_per-sample.png",
  width = 7, height = 5
  )
dimplot(
  obj       = integrated_rpca,
  group_by  = "integrated_annotation",
  split_by  = "orig.ident",
  pallette  = get_color_pallette("per_celltype"),
  file_path = "./Root/UMAP plots/split_integrated_RPCA_UMAP_per_celltype.png"
  )
dimplot(
  obj       = integrated_rpca,
  group_by  = "integrated_annotation",
  split_by  = "orig.ident",
  pallette  = get_color_pallette("per_cluster"),
  file_path = "./Root/UMAP plots/split_integrated_RPCA_UMAP_per_cluster.png"
  )