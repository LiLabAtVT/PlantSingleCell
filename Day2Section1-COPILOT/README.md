# load conda module
module load Miniconda3/24.7.1-0

# activate conda environment
source activate /projects/songli_lab/PlantSingleCell2025/Day_2/scKB_copilot



# inside R - load R environment
Sys.setenv(R_LIBS_USER="/projects/songli_lab/PlantSingleCell2025/Day_2/scKB_copilot_R")
.libPaths(Sys.getenv("R_LIBS_USER"))
.libPaths()  # Confirm the change