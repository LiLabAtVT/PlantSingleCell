# Create and go to your own folder
mkdir /projects/songli_lab/PlantSingleCell2025/Day_2/ParticipantFolder/{your_name} \
cd /projects/songli_lab/PlantSingleCell2025/Day_2/ParticipantFolder/{your_name}

# Download the data, code and supplementary files and prepare them for use
git clone https://github.com/ohlerlab/scKB.git 

cd ./scKB \
unzip 10xv2_whitelist.txt.zip # cell barcode white list for 10x chemistry 2 \
unzip 10xv3_whitelist.txt.zip # cell barcode white list for 10x chemistry 3 \
gunzip Arabidopsis_thaliana.TAIR10.43.gtf.gz \
chmod 777 scKB

# Load conda module
module load Miniconda3/24.7.1-0

# Activate conda environment
source activate /projects/songli_lab/PlantSingleCell2025/Day_2/scKB_copilot

# Inside R - load R environment
R \
Sys.setenv(R_LIBS_USER="/projects/songli_lab/PlantSingleCell2025/Day_2/scKB_copilot_R") \
.libPaths(Sys.getenv("R_LIBS_USER")) \
.libPaths()  # Confirm the change

# Prepare the genome for alignment with kallisto and bustools in R
library(BUSpaRse) \
library(BSgenome.Athaliana.TAIR.TAIR9) # Load the Arabidopsis genome \
get_velocity_files(X = "./Arabidopsis_thaliana.TAIR10.43.gtf", L = 91, Genome = BSgenome.Athaliana.TAIR.TAIR9, out_path = "./", isoform_action = "separate", chrs_only=FALSE, style="Ensembl") # gtf is the genome annotation file, and 91 is the R2 length of the toy data (10x chemistry v3) \
system("kallisto index -i ./cDNA_introns_10xv3.idx ./cDNA_introns.fa") # Index the intron file \
quit()

# Align the single cell raw reads to genome with kallisto and bustools
./scKB -f ./toy_data -i ./cDNA_introns_10xv3.idx -d ./ -s 10xv3 -t 16 -w ./10xv3_whitelist.txt -n ./col0_toy

# Run COPILOT
R \
Sys.setenv(R_LIBS_USER="/projects/songli_lab/PlantSingleCell2025/Day_2/scKB_copilot_R") \
.libPaths(Sys.getenv("R_LIBS_USER")) \
.libPaths() # Confirm the change \
library(COPILOT) \
copilot(sample.name = "col0_toy", species.name = "Arabidopsis thaliana", transcriptome.name = "TAIR10", sample.stats = NULL, mt.pattern = "ATMG", mt.threshold = 5, cp.pattern = "ATCG", remove.doublet = FALSE, do.seurat = FALSE, do.annotation= FALSE, unwanted.genes = NULL, filtering.ratio = 1, min.UMI.low.quality = 1, min.UMI.high.quality = 3) \
quit() ## Check out the col0_toy for results!


## Please look up the parameters and more documentations at https://github.com/Hsu-Che-Wei/COPILOT

## Here is the full codes to create the scKB_copilot environment from scratch (The one described in the STARProtocol published back in 2022 is already outdated):

conda create -n scKB_copilot -c conda-forge -c bioconda r-base=4.1.3 kallisto=0.48.0 bustools=0.41.0 bioconductor-busparse=1.8.0 bioconductor-bsgenome=1.62.0 bioconductor-dropletutils=1.14.2 r-devtools=2.4.3 r-rjson=0.2.21 r-r2html=2.3.2 python=3.8 jupyterlab pyscaffold -y \
conda activate scKB_copilot 

conda install conda-forge::r-ape -y \
conda install conda-forge::r-reticulate -y \
conda install conda-forge::umap-learn -y \
conda install conda-forge::leidenalg -y \
conda install anaconda::pandas -y \
conda install jupyterlab -y 

R #Do not update old packages!! \
install.packages("BiocManager") \
install.packages(c('repr', 'IRdisplay', 'evaluate', 'crayon', 'pbdZMQ', 'uuid', 'digest', 'IRkernel')) \
BiocManager::install("BSgenome.Athaliana.TAIR.TAIR9") \
install.packages(c('R2HTML','remotes','cluster', 'cowplot', 'fitdistrplus', 'future', 'future.apply', 'ggrepel', 'ggridges', 'ica', 'igraph', 'irlba', 'leiden', 'lmtest', 'patchwork', 'pbapply', 'plotly', 'RANN', 'RcppAnnoy', 'reticulate', 'tsne', 'Rtsne', 'sctransform', 'uwot', 'RcppEigen','rsvd')) \
remotes::install_github('https://github.com/ekernf01/DoubletFinder', force = T) \
install.packages("/projects/songli_lab/PlantSingleCell2025/Day_2/ParticipantFolder/CheWei/Seurat_3.1.5.tar.gz") \
devtools::install_github("ohlerlab/COPILOT") \
IRkernel::installspec()




