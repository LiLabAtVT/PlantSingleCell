#!/bin/bash
#SBATCH --job-name=FilterConvert     # Job name
#SBATCH --nodes=1                                # Number of nodes
#SBATCH --ntasks=1                               # Number of tasks (MPI processes)
#SBATCH --cpus-per-task=4                        # Number of CPU cores per task
#SBATCH --mem=16G                                # Memory per node (in GB)
#SBATCH --time=20:00:00                          # Time limit (hh:mm:ss)
#SBATCH --mail-user=arazan@vt.edu                # Email address for job notifications
#SBATCH --mail-type=ALL                          # Send email at beginning and end of job
#SBATCH --output=filter_convert_Matthew_10X_%j.out           # Standard output log
#SBATCH --account=introtogds                     # Replace with your valid account

# Define paths to custom tools
BAMTOFASTQ="/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/Convert_BAM_to_FASTQ/bamtofastq_linux"
SUBSET_BAM="/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/Convert_BAM_to_FASTQ/subset-bam_linux"

# Define input files and output directories with unique timestamps
input_bam="possorted_genome_Pos24hpi_1_bam.bam"        # Main BAM file
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
