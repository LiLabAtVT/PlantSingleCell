# 🌱 PlantSingleCell 2025 – Day 1 Session 2 Setup Guide

Welcome to **Session 2** of the PlantSingleCell 2025 Workshop! This guide walks you through setting up your environment and starting your analysis using RStudio on VT ARC’s OOD platform.

---

## 🎯 Objective

The goal of this session is to help you set up a reproducible and high-performance computing environment for integrating single-cell plant data using R. You will learn how to:
- Access and clone the workshop repository
- Launch and configure RStudio on the Virginia Tech ARC cluster
- Set up a dedicated R environment and working directory
- Begin analysis using scripts provided in the `Day1Section2 (Integration)` folder


---

## 📁 Step 1: Get the Scripts

1. **Create your working directory**:
   ```bash
   mkdir -p /projects/songli_lab/PlantSingleCell2025/Day_1/Session_2/[dir_name]
   cd /projects/songli_lab/PlantSingleCell2025/Day_1/Session_2/[dir_name]
   ```

2. **Clone the repository**:
   ```bash
   git clone https://github.com/LiLabAtVT/PlantSingleCell.git
   ```

3. **Navigate into the Integration directory**:
   ```bash
   cd PlantSingleCell/Day1Section2 (Integration)
   ```

---

## 💻 Step 2: Launch RStudio on VT ARC

Go to 👉 [https://ood.arc.vt.edu/pun/sys/dashboard](https://ood.arc.vt.edu/pun/sys/dashboard) and click on RStudio.

### 📌 RStudio Launch Settings:
- **Cluster**: OWL - Water-cooled AMD CPU  
- **Partition**: `normal_q`  
- **Account**: `introtogds`  
- **Run Time**: `2 hours`  
- **Version**: RStudio 2024.12.0 / R 4.4.2  
- **Cores**: `12`  

### ⚙️ Advanced Options:
- **Override Memory (GB)**: `256`  
- **Extra Modules**:
  ```
  GSL/2.7-GCC-13.2.0
  HDF5/1.14.3-gompi-2023b
  ```

👉 Click **Launch**, and once the status changes to **Running**, hit **“Connect to RStudio Server”**.

---

## 🧪 Step 3: Prepare the R Environment

In the **Console tab** of RStudio, run the following:

1. **Set your R library path**:
   ```r
   Sys.setenv(R_LIBS_USER="/projects/songli_lab/PlantSingleCell2025/Day_1/Session_2/env/")
   .libPaths(Sys.getenv("R_LIBS_USER"))
   .libPaths()  # Confirm the change
   ```

2. **Set the working directory**:
   ```r
   setwd("/projects/songli_lab/PlantSingleCell2025/Day_1/Session_2/[dir_name]/PlantSingleCell/Day1Section2 (Integration)")
   ```

---

## 🚀 Step 4: Start Your Analysis

You're all set! 🎉 Start running the commands provided in the scripts and begin exploring your data.

Happy analyzing! 🍀
