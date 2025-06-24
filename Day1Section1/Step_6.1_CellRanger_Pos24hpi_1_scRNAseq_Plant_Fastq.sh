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

