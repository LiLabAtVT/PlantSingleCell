# Plant Single Cell Analysis Workshop and Mini-hackathon
This is the github repository for a summer training camp for plant single cell analysis at Virginia Tech


## Day 1 Materials
### Day 1 section 1 (Razan) Read mapping with ARC, and load data into Seurat.
[View detailed setup instructions](./Day1Section1/README.md)

### Day 1 section 2 (Maryam) Integration across different runs, etc.
[View detailed setup instructions](./Day1Section2-Integration/README.md)

## Day 2 Materials
### Day 2 section 1 (CoPilot method)
Dr. Wei-Che Hsu will present COPILOT[https://github.com/Hsu-Che-Wei/COPILOT]

### Day 2 section 2 (Michael Passalacqua)
Mr. Michael Passalacqua from Cold Spring Harbor lab will present his co-expression proxy method[https://github.com/gillislab/Coexpression_Proxies]

### Day 2 section 3 
Ms. Tran Chau will present Othro Marker Gene groups method.
To quickly install the required packages, you need to run this
```bash
Sys.setenv(R_LIBS_USER="/projects/songli_lab/PlantSingleCell2025/Day_2/env/")
.libPaths(Sys.getenv("R_LIBS_USER"))
.libPaths() 
```
Use OMG for cell type annotation, https://github.com/LiLabAtVT/OrthoMarkerGeneGroups

### Day 2 section 5 (Ms. Tran Chau will present the SATURN method walk-through)
