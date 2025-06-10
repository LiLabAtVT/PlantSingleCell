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
CELLRANGER_PATH="/projects/intro2gds/Razan_2024/scRNA_Seq_Arab/Refdata_Arab/cellranger-8.0.1"  # Update this if needed

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
