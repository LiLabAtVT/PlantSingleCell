# Conda Environment Setup on VT ARC HPC

This guide walks you through setting up and using a **Conda Python environment** on the Virginia Tech ARC system using the provided Miniconda module.

---

## Step 1: Load Miniconda Module

ARC provides Miniconda via modules. Load it like this:

```bash
module load Miniconda3/24.7.1-0
```

---

## Step 2: Create and Activate a Conda Environment

Create a new environment named `myenv` with Python 3.9:

```bash
conda create --name myenv python=3.9
```

Activate it:

```bash
conda activate myenv
```

Install commonly used packages:

```bash
conda install numpy pandas matplotlib
```

---

## Step 3 (Optional): Create from environment.yml

If you already have an environment file:

```bash
conda env create -f environment.yml
conda activate myenv
```

---

## Step 4: Save or Remove Environment

Export your environment:

```bash
conda env export > environment.yml
```

Remove an environment:

```bash
conda remove --name myenv --all
```

---

## Tips

- Conda environments are stored in your `$HOME`, so they persist across sessions.
- Use named environments instead of modifying the base environment.
- Deactivate any active environment with:

```bash
conda deactivate
```

---

## TL;DR Quick Start

```bash
module load Miniconda3/24.7.1-0
conda create -n myenv python=3.9
conda activate myenv
conda install numpy pandas
```

---

## Resources

- Conda Docs: https://docs.conda.io
- VT ARC HPC: https://arc.vt.edu
