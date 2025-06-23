# 🌱 Plant Single Cell Analysis Workshop & Mini-Hackathon
**📍 Virginia Tech | Summer Training Camp 🧡🦃**

Welcome to the official repository for the **Plant Single Cell Analysis** summer camp! This 4-day workshop is designed to equip participants with cutting-edge tools, workflows, and insights into plant single-cell transcriptomics.

---
### 🚀 Quick Start Instructions

- [📦 Create Conda Environment](./Create_Conda_env)
- [📦 Launch Apptainer Demo](./Launch_Apptainer)
- [📦 Launch RStudio Guide](./Launch_RStudio)

> Click on each folder to view the step-by-step setup instructions and example files.
---

## 📅 Day 1: Foundations of Single-Cell Analysis

### 🧬 Section 1: Read Mapping & Seurat Integration
**Presenter: Razan Alajoleen**

🧾 Topics Covered:
- Read mapping using ARC
- Loading single-cell data into Seurat

🔗 [Click here for detailed setup instructions](./Day1Section1/README.md)

---

### 🧬 Section 2: Cross-Sample Integration with Seurat
**Presenter: Maryam Haghani**

🧾 Topics Covered:
- Sample integration
- UMAP visualization: before vs. after integration
- Gene co-expression of marker genes in *Arabidopsis* root

🔗 [Click here for detailed setup instructions](./Day1Section2-Integration/README.md)

---

## 📅 Day 2: Advanced Methods and Tools

### 🧠 Section 1: COPILOT – Cell Type Classification with Deep Learning
**Presenter: Dr. Wei-Che Hsu**

<img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/github/github-original.svg" width="20"/> [Explore COPILOT on GitHub](https://github.com/Hsu-Che-Wei/COPILOT)

---

### 🔬 Section 2: Co-expression Proxy Method
**Presenter: Mr. Michael Passalacqua (Cold Spring Harbor Lab)**

<img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/github/github-original.svg" width="20"/> [Explore Coexpression_Proxies on GitHub](https://github.com/gillislab/Coexpression_Proxies)

---

### 🌿 Section 3: Orthologous Marker Gene Groups (OMG)
**Presenter: Tran Chau**

🧾 Topics Covered:
- Cross-species marker gene identification
- Cell type annotation using orthologous gene sets

🔧 To use required packages:
```r
Sys.setenv(R_LIBS_USER="/projects/songli_lab/PlantSingleCell2025/Day_2/env/")
.libPaths(Sys.getenv("R_LIBS_USER"))
.libPaths()
```

<img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/github/github-original.svg" width="20"/> [Explore OMG on GitHub](https://github.com/LiLabAtVT/OrthoMarkerGeneGroups)


---

### 🔁 Section 5: SATURN Method Walk-through
**Presenter: Tran Chau**

🧾 A step-by-step overview of the SATURN method for cell annotation.

<img src="https://cdn.jsdelivr.net/gh/devicons/devicon/icons/github/github-original.svg" width="20"/> [Explore SATURN on GitHub](https://github.com/snap-stanford/SATURN)

🔗 [Click here for detailed setup instructions](./Saturn_Plant/README.md)


---

We hope this workshop helps you discover the exciting frontier of plant single-cell genomics! 🌾

---

## 👥 Organizers

### 🔬 Principal Investigator
- **[Song Li](https://github.com/songliVT)** – Associate Professor, Virginia Tech  
  Lead organizer and PI for the Plant Single Cell Analysis Workshop

### 🧑‍💻 Core Organizing Team
- **[Razan Alajoleen (Postdoctoral Research Associate)](https://github.com/RazanAj)** 
- **[Tran N. Chau (Ph.D. Student)](https://github.com/ct-tranchau)**  
- **[Maryam Haghani (Ph.D. Student)](https://github.com/Maryam-Haghani)** 
- **[Sajib Acharjee Dip (Ph.D. Student)](https://github.com/Sajib-006)**

> Special thanks to all contributors for making the Summer Training Camp a success! 🎉

---

