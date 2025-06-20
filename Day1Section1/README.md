# 📘 Day 1 – Session 1: Read Mapping for Plant Single-Cell Analysis  
**Data: Pathogen-Infected Arabidopsis Roots** (snRNA seq Data)

# 🧬 Objective  
To perform read alignment and gene counting for a single-cell RNA-seq dataset from pathogen-infected plant roots using 10X Genomics data. The example uses Arabidopsis thaliana as the host.

# Quick Start

1. **Create your working directory**:
   ```bash
   mkdir -p /projects/songli_lab/PlantSingleCell2025/Day_1/ParticipantFolder/[dir_name]
   cd /projects/songli_lab/PlantSingleCell2025/Day_1/ParticipantFolder/[dir_name]
   ```

2. **Clone the repository**:
   ```bash
   git clone https://github.com/LiLabAtVT/PlantSingleCell.git
   ```

3. **Navigate into the Section 1 directory**:
   ```bash
   cd PlantSingleCell/Day1Section1
   ```
4. **Submit your SLURM script**:
   ```bash
   sbatch Cellranger_mkref_Multiplespecies.sh
   ```

***Steps Covered***

### Build a Combined Reference

1. **Download FASTA and GTF files**  
   - Arabidopsis genome and annotation from [Ensembl Plants](https://plants.ensembl.org/Arabidopsis_thaliana/Info/Index)
   - Pathogen genome (Phytophthora capsici) from NCBI (https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_016618375.1/)

2. **Build Cell Ranger-Compatible Reference**  
- Use makeref command to create a reference (https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-3p-references#multiple-species-4f40e4)

3. **Download Cell Ranger**
 https://www.10xgenomics.com/support/software/cell-ranger/downloads#download-links
  
```
#!/bin/bash
#SBATCH --job-name=cellranger_mkref_multi      # Job name
#SBATCH --nodes=1                              # Number of nodes
#SBATCH --ntasks=1                             # Number of tasks (MPI processes)
#SBATCH --cpus-per-task=4                      # Number of CPU cores per task
#SBATCH --mem=50G                              # Memory per node (in GB)
#SBATCH --time=12:00:00                        # Time limit (hh:mm:ss)
#SBATCH --mail-user=arazan@vt.edu              # Email address for job notifications
#SBATCH --mail-type=ALL                        # Send email at beginning and end of job
#SBATCH --output=Cellranger_mkref_%j.out       # Standard output log
#SBATCH --account=introtogds                   # Replace with your valid account

echo "Starting cellranger mkref for multiple species"
date

# Full path to the cellranger executable
CELLRANGER_PATH="/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/cellranger-9.0.1"  # Update this if needed

# Define paths for multiple species
GENOME1_NAME="Arabidopsis_thaliana"
FASTA_FILE1="/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Ref_Genome10X_Multiple_Species_Genome/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa"
GTF_FILE1="/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Ref_Genome10X_Multiple_Species_Genome/Arabidopsis_thaliana.TAIR10.59.gtf.filtered.gtf"

GENOME2_NAME="Pcap"
FASTA_FILE2="/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Ref_Genome10X_Multiple_Species_Genome/GCA_016618375.1_Pcap_4.1_genomic_cleaned.fasta"
GTF_FILE2="/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Ref_Genome10X_Multiple_Species_Genome/GCA_016618375.1_Pcap_4.1_genomic_Filtered.gtf"

# Check if files exist for both species
for file in "$FASTA_FILE1" "$GTF_FILE1" "$FASTA_FILE2" "$GTF_FILE2"; do
    if [ ! -f "$file" ]; then
        echo "Error: File $file does not exist."
        exit 1
    fi
done

# Run cellranger mkref for multiple species
$CELLRANGER_PATH/cellranger mkref \
    --genome=$GENOME1_NAME \
    --fasta=$FASTA_FILE1 \
    --genes=$GTF_FILE1 \
    --genome=$GENOME2_NAME \
    --fasta=$FASTA_FILE2 \
    --genes=$GTF_FILE2

# Check if cellranger mkref succeeded
if [ $? -ne 0 ]; then
    echo "Error: cellranger mkref command failed."
    exit 1
fi

echo "cellranger mkref for multiple species completed successfully"
date

exit 0
```
### Run `cellranger count` to generate BAM and Barcode CSV Files

Use `cellranger count` to generate:
- A BAM file containing aligned reads
- A CSV file with barcodes assigned to each species (plant and pathogen)

This step uses the **Build a Combined Reference** created in the previous step.
```
#!/bin/bash
#SBATCH --job-name=Cellranger_Analysis_Pos24hpi_1  # Job name
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of tasks (MPI processes)
#SBATCH --cpus-per-task=64               # Number of CPU cores per task
#SBATCH --mem=128G                       # Memory per node (in GB)
#SBATCH --time=24:00:00                 # Time limit (hh:mm:ss)
#SBATCH --mail-type=ALL                      # Send email at beginning and end of job
#SBATCH --mail-user=arazan@vt.edu       # Email address for job notifications
#SBATCH --output=Cellranger_Pos24hpi_Multi_Ref.out      # Standard output log
#SBATCH --account=introtogds            # Replace with your valid account

echo "Starting"

# Run Cell Ranger count command
/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/cellranger-9.0.1/cellranger count \
    --id=Plant_Pathogen_Mapping_Pos24hpi_1_scRNAseq_Multi_Ref \
    --transcriptome=/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Ref_Genome10X_Multiple_Species_Genome/Arabidopsis_thaliana_and_Pcap \
    --fastqs=/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/Pos24hpi_1 \
    --sample=Pos24hpi_1_CKDL240032444-1A_22GH3VLT4 \
    --force-cells=8000 \
    --create-bam=true


echo "Finished"
date

exit;
```
```
#!/bin/bash
#SBATCH --job-name=Cellranger_Analysis_Neg24hpi_1  # Job name
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of tasks (MPI processes)
#SBATCH --cpus-per-task=64               # Number of CPU cores per task
#SBATCH --mem=128G                       # Memory per node (in GB)
#SBATCH --time=24:00:00                 # Time limit (hh:mm:ss)
#SBATCH --mail-type=ALL                      # Send email at beginning and end of job
#SBATCH --mail-user=arazan@vt.edu       # Email address for job notifications
#SBATCH --output=Cellranger_Neg24hpi_Multi_Ref.out      # Standard output log
#SBATCH --account=introtogds            # Replace with your valid account

echo "Starting"

# Run Cell Ranger count command
/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/cellranger-9.0.1/cellranger count \
    --id=PLant_Pathogen_Mapping_Neg24hpi_1_scRNAseq \
    --transcriptome=/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Ref_Genome10X_Multiple_Species_Genome/Arabidopsis_thaliana_and_Pcap \
    --fastqs=/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/Neg24hpi_1 \
    --sample=Neg24hpi_1_CKDL240032442-1A_22GH3VLT4 \
    --force-cells=8000 \
    --create-bam=true
  

echo "Finished"
date

exit;

```


### Identify Barcodes 
Make sure the gem_classification.csv files for both Pos24 and Neg24 samples are present in your working directory. These files can be found in the respective analysis folders generated by the Cell Ranger output.

Run the following command: 
```
awk -F',' 'NR > 1 && $4 == "Arabidopsis_thaliana" { print $1 }' gem_classification_Pos24hpi_1.csv > arabidopsis_barcodes.txt 
awk -F',' 'NR > 1 && $4 == "Pcap" { print $1 }' gem_classification_Pos24hpi_1.csv > pcap_barcodes.txt 
awk -F, '$4 == "Multiplet" { print $1 }' gem_classification_Pos24hpi_1.csv > multiplet_barcodes.txt
```
### Generate BAM Files- Convert BAM to FASTQ
Make sure the possorted_genome_bam.bam and the possorted_genome-bam.bam.bai files are present in your working directory. These files are in the OUT folder generated by the Cell Ranger output.

*Tools needed:
1. bamtofastq_linux
2. subset-bam_linux
Function: Splits BAM by barcode and converts to FASTQ

https://github.com/10XGenomics/subset-bam/releases 

https://github.com/10XGenomics/bamtofastq/releases 

3. Samtools:
   https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2

```
#!/bin/bash
#SBATCH --job-name=FilterConvert     # Job name
#SBATCH --nodes=1                                # Number of nodes
#SBATCH --ntasks=1                               # Number of tasks (MPI processes)
#SBATCH --cpus-per-task=4                        # Number of CPU cores per task
#SBATCH --mem=16G                                # Memory per node (in GB)
#SBATCH --time=20:00:00                          # Time limit (hh:mm:ss)
#SBATCH --mail-user=arazan@vt.edu                # Email address for job notifications
#SBATCH --mail-type=ALL                          # Send email at beginning and end of job
#SBATCH --output=filter_convert_%j.out           # Standard output log
#SBATCH --account=introtogds                     # Replace with your valid account

# Define paths to custom tools
BAMTOFASTQ="/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/Convert_BAM_to_FASTQ/bamtofastq_linux"
SUBSET_BAM="/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/Convert_BAM_to_FASTQ/subset-bam_linux"

# Define input files and output directories with unique timestamps
input_bam="/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/Convert_BAM_to_FASTQ/possorted_genome_Pos24hpi_1_bam.bam"  # Main BAM file
plant_barcodes="arabidopsis_barcodes.txt"              # Plant barcodes file
pathogen_barcodes="pcap_barcodes.txt"                  # Pathogen barcodes file
plant_bam="plant_reads.bam"                            # Filtered plant BAM
pathogen_bam="pathogen_reads.bam"                      # Filtered pathogen BAM
timestamp=$(date +"%Y%m%d_%H%M%S")                     # Unique timestamp for output
plant_fastq_dir="./plant_fastq_output_$timestamp"      # Unique directory for plant FASTQ files
pathogen_fastq_dir="./pathogen_fastq_output_$timestamp" # Unique directory for pathogen FASTQ files

# Step 1: Index the main BAM file if not already indexed
if [ ! -f "${input_bam}.bai" ]; then
    echo "Indexing BAM file..."
    samtools index "$input_bam"
fi

# Step 2: Remove existing BAM files if they exist
if [ -f "$plant_bam" ]; then
    echo "Removing existing plant BAM file..."
    rm -f "$plant_bam"
fi

if [ -f "$pathogen_bam" ]; then
    echo "Removing existing pathogen BAM file..."
    rm -f "$pathogen_bam"
fi

# Step 3: Filter BAM file for plant and pathogen reads using subset-bam
echo "Filtering BAM file for plant-specific reads..."
$SUBSET_BAM --bam "$input_bam" --cell-barcodes "$plant_barcodes" --out-bam "$plant_bam"

echo "Filtering BAM file for pathogen-specific reads..."
$SUBSET_BAM --bam "$input_bam" --cell-barcodes "$pathogen_barcodes" --out-bam "$pathogen_bam"

# Step 4: Convert filtered BAM files to FASTQ using 10x's bamtofastq
echo "Converting plant BAM to FASTQ..."
#mkdir -p "$plant_fastq_dir"
$BAMTOFASTQ "$plant_bam" "$plant_fastq_dir"

echo "Converting pathogen BAM to FASTQ..."
#mkdir -p "$pathogen_fastq_dir"
$BAMTOFASTQ "$pathogen_bam" "$pathogen_fastq_dir"

echo "Process complete! FASTQ files for plant and pathogen are saved in their respective directories."
```
### Process FASTQ Files_Run Cell Ranger Count 

```
#!/bin/bash
#SBATCH --job-name=Cellranger_Analysis_Pos24hpi_1  # Job name
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of tasks (MPI processes)
#SBATCH --cpus-per-task=64               # Number of CPU cores per task
#SBATCH --mem=128G                       # Memory per node (in GB)
#SBATCH --time=24:00:00                 # Time limit (hh:mm:ss)
#SBATCH --mail-type=ALL                      # Send email at beginning and end of job
#SBATCH --mail-user=arazan@vt.edu       # Email address for job notifications
#SBATCH --output=Cellranger_Pos24hpi_1_ATH_Read_Only.out      # Standard output log
#SBATCH --account=introtogds            # Replace with your valid account

echo "Starting"

# Run Cell Ranger count command
/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/cellranger-9.0.1/cellranger count \
    --id=Mapping_Pos24hpi_1_scRNAseq_ATH_Read_Only \
    --transcriptome=/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Ref_Genome10X_Plant/ATH_genome \
    --fastqs=/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/Convert_BAM_to_FASTQ/plant_fastq_output_20250613_114515/Plant_Pathogen_Mapping_Pos24hpi_1_scRNAseq_Multi_Ref_0_1_22GH3VLT4 \
    --sample=Pos24hpi1 \
    --force-cells=8000 \
    --create-bam=true

echo "Finished"
date
exit;

```
### Load Filtered Raw Data into RStudio and Perform Quality Control (QC)

Load the filtered raw gene expression data into RStudio to inspect quality control (QC) metrics using the script below.

```
library(Seurat)
library(sctransform)

# Load raw 10X count matrices
Neg1 <- Read10X(data.dir = "/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/CellRanger_Output_Mapped_Multi_Ref/PLant_Pathogen_Mapping_Neg24hpi_1_scRNAseq_Multi_Ref/outs/filtered_feature_bc_matrix")
Pos1 <- Read10X(data.dir = "/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/Convert_BAM_to_FASTQ/Mapping_Pos24hpi_1_scRNAseq_ATH_Read_Only/outs/filtered_feature_bc_matrix")

rownames(Neg1) <- gsub("Arabidopsis_thaliana_gene:", "", rownames(Neg1))
rownames(Pos1) <- gsub("gene:", "", rownames(Pos1))


process_seurat_object <- function(counts_data, project_name,
                                  feature_upper, nfeatures,
                                  pca_dims, clustering_resolution) {
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts_data, project = project_name, 
                                   min.cells = 3, min.features = 200)
  
  # Compute mitochondrial percentage for Arabidopsis
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^ATMG")

  # Arabidopsis: no usable percent.mt, so set it to 0
  seurat_obj[["percent.mt"]] <- 0
  
  # Filter based on quality control metrics
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 100 & 
                         nFeature_RNA < feature_upper & 
                         percent.mt < 5)
   # QC violin plots
  pdf(paste0("./", project_name, "_QC_violin.pdf"), width = 10, height = 5)
  print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  dev.off()  
  
  # Normalize data
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  # Identify variable features
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = nfeatures)
  # Scale and run PCA
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  # Cluster and visualize
  seurat_obj <- FindNeighbors(seurat_obj, dims = pca_dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = clustering_resolution)
  seurat_obj <- RunUMAP(seurat_obj, dims = pca_dims)
  
  # Save object
  saveRDS(seurat_obj, paste0("./", project_name, ".rds"))
  
  return(seurat_obj)
}


Neg24hpi_1 <- process_seurat_object(Neg1, project_name = "Neg_hpi1", 
                                    feature_upper = 6000, nfeatures = 2000,
                                    pca_dims = 1:20, clustering_resolution = 0.22)

Pos24hpi_1 <- process_seurat_object(Pos1, project_name = "Pos_hpi1", 
                                    feature_upper = 3500, nfeatures = 2000,
                                    pca_dims = 1:20, clustering_resolution = 0.22)




crossSamples <- SplitObject(merge(Neg24hpi_1, y= c(Pos24hpi_1)) , split.by = "orig.ident")
for (i in names(crossSamples)) {
  crossSamples[[i]] <- SCTransform(crossSamples[[i]], verbose = FALSE)
}
features <- SelectIntegrationFeatures(object.list = crossSamples, nfeatures = 10000)
print(length(features))
crossSamples <- PrepSCTIntegration(object.list = crossSamples, anchor.features = features)
crossSamples.anchors <- FindIntegrationAnchors(object.list = crossSamples, normalization.method = "SCT", anchor.features = features)
crossSamples.combine <- IntegrateData(anchorset = crossSamples.anchors, normalization.method = "SCT")
crossSamples.combine <- RunPCA(crossSamples.combine, verbose = FALSE)
crossSamples.combine <- RunUMAP(crossSamples.combine, reduction = "pca", dims = 1:30)
crossSamples.combine <- FindNeighbors(crossSamples.combine, reduction = "pca", dims = 1:30)
crossSamples.combine <- FindClusters(crossSamples.combine, resolution = 0.3) #org =0.3

saveRDS(crossSamples.combine, "./Patho_Nonpatho_integration.Rds")

pdf("./Patho_Nonpatho.pdf", width = 8, height = 5)
DimPlot(crossSamples.combine, reduction = "umap", group.by = "orig.ident", cols = c("lightblue","blue", "#FFC073","#FF8C00"))
dev.off()

pdf("./Patho_Nonpatho_sepColor.pdf", width = 8, height = 5)
DimPlot(crossSamples.combine, reduction = "umap", split.by = "orig.ident", group.by = "orig.ident", cols = c("lightblue","blue", "#FFC073","#FF8C00"))
dev.off()

pdf("./Patho_Nonpatho_lable.pdf", width = 8, height = 5)
DimPlot(crossSamples.combine, reduction = "umap", label = TRUE)
dev.off()


pdf("./Patho_Nonpatho_sep.pdf", width = 8, height = 5)
DimPlot(crossSamples.combine, reduction = "umap", split.by = "orig.ident")
dev.off()

```
