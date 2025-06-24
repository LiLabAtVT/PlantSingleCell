#!/bin/bash
#SBATCH --job-name=cellranger_mkref           # Job name
#SBATCH --nodes=1                             # Number of nodes
#SBATCH --ntasks=1                            # Number of tasks (MPI processes)
#SBATCH --cpus-per-task=4                     # Number of CPU cores per task
#SBATCH --mem=50G                             # Memory per node (in GB)
#SBATCH --time=12:00:00                       # Time limit (hh:mm:ss)
#SBATCH --mail-user=arazan@vt.edu             # Email address for job notifications
#SBATCH --mail-type=ALL                       # Send email at beginning and end of job
#SBATCH --output=cellranger_mkref_%j.out      # Standard output log
#SBATCH --account=introtogds                  # Replace with your valid account

echo "Starting cellranger mkref"
date

# Full path to the cellranger executable
CELLRANGER_PATH="/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Data/cellranger-9.0.1"  # Replace with the actual path

# Set paths for the genome and GTF files
GENOME_NAME="ATH_genome"
FASTA_FILE="/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Ref_Genome10X_Plant/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa"
GTF_FILE="/projects/songli_lab/PlantSingleCell2025/Day_1/Session_1/Ref_Genome10X_Plant/Arabidopsis_thaliana.TAIR10.59.gtf.filtered.gtf"

# Check if the FASTA file exists
if [ ! -f "$FASTA_FILE" ]; then
    echo "Error: FASTA file $FASTA_FILE does not exist."
    exit 1
fi

# Check if the GTF file exists
if [ ! -f "$GTF_FILE" ]; then
    echo "Error: GTF file $GTF_FILE does not exist."
    exit 1
fi

# Run cellranger mkref
$CELLRANGER_PATH/cellranger mkref \
    --genome=$GENOME_NAME \
    --fasta=$FASTA_FILE \
    --genes=$GTF_FILE 

#if [ $? -ne 0 ]; then
   # echo "Error: cellranger mkref command failed."
    #echo "Saving diagnostic information..."
    #if [ -f "mkref_${GENOME_NAME}/mkref_${GENOME_NAME}.mri.tgz" ]; then
       # $CELLRANGER_PATH/cellranger upload arazan@tinkercliffs1.arc.vt.edu "mkref_${GENOME_NAME}/mkref_${GENOME_NAME}.mri.tgz"
    #else
        #echo "Error: Diagnostic file not found."
    #fi
    #exit 1
#fi

echo "cellranger mkref completed successfully"
date

exit 0
