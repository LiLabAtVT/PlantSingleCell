# load conda module
module load Miniconda3/24.7.1-0

# activate conda environment
source activate /projects/songli_lab/PlantSingleCell2025/Day_2/scKB_copilot

# inside R - load R environment
Sys.setenv(R_LIBS_USER="/projects/songli_lab/PlantSingleCell2025/Day_2/scKB_copilot_R") \
.libPaths(Sys.getenv("R_LIBS_USER")) \
.libPaths()  # Confirm the change

# Create and go to your own folder
mkdir /projects/songli_lab/PlantSingleCell2025/Day_2/ParticipantFolder/{your_name} \
cd /projects/songli_lab/PlantSingleCell2025/Day_2/ParticipantFolder/{your_name}

# Download the data, code and supplementary files and prepare them for use
git clone https://github.com/ohlerlab/scKB.git \

cd ./scKB \
unzip 10xv2_whitelist.txt.zip \
unzip 10xv3_whitelist.txt.zip \
gunzip Arabidopsis_thaliana.TAIR10.43.gtf.gz \
chmod 777 scKB

# Prepare for the genome for alignment with kallisto and bustools in R
R \ 
library(BUSpaRse) \
library(BSgenome.Athaliana.TAIR.TAIR9) # Load the Arabidopsis genome \
get_velocity_files(X = "./Arabidopsis_thaliana.TAIR10.43.gtf", L = 91, Genome = BSgenome.Athaliana.TAIR.TAIR9, out_path = "./", isoform_action = "separate", chrs_only=FALSE, style="Ensembl") # gtf is the genome annotation file, and 91 is the R2 length of the toy data (10x chemistry v3) \
system("kallisto index -i ./cDNA_introns_10xv3.idx ./cDNA_introns.fa") # Index the intron file \
quit()

# Align the single cell raw reads to genome with kallisto and bustools
./scKB -f ./toy_data -i ./cDNA_introns_10xv3.idx -d ./ -s 10xv3 -t 16 -w ./10xv3_whitelist.txt -n ./col0_toy

# Run COPILOT
R \
library(COPILOT) \
copilot(sample.name = "col0_toy", species.name = "Arabidopsis thaliana", transcriptome.name = "TAIR10", sample.stats = NULL, mt.pattern = "ATMG", mt.threshold = 5, cp.pattern = "ATCG", remove.doublet = FALSE, do.seurat = FALSE, do.annotation= FALSE, unwanted.genes = NULL, filtering.ratio = 1, min.UMI.low.quality = 1, min.UMI.high.quality = 3)

## Please look up the paramenters at https://github.com/Hsu-Che-Wei/COPILOT





