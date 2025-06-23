# Using Apptainer on VT ARC HPC Systems

This guide helps you set up and run containers using **Apptainer** (formerly Singularity) on Virginia Tech’s ARC clusters via terminal or Open OnDemand.

---

## Step 1: Load Apptainer

Apptainer is provided via the module system.

```bash
module load apptainer/1.4.0
```

Verify the installation:

```bash
apptainer --version
```

Expected output:

```
apptainer version 1.4.0
```

---

## Step 2: Run a Simple Test

Run a quick test using a lightweight container:

```bash
apptainer run docker://alpine echo "Hello from inside a container!"
```

Expected output:

```
Hello from inside a container!
```

---

## Step 3: Shell into a Container

You can get an interactive shell in a container like this:

```bash
apptainer shell docker://ubuntu:20.04
```

Try running some commands inside the container:

```bash
ls /
exit
```

---

## Step 4 (Optional): Build a `.sif` Container Image

If building is allowed on your node:

```bash
apptainer build ubuntu_20.04.sif docker://ubuntu:20.04
```

Then run offline using:

```bash
apptainer run ubuntu_20.04.sif echo "Offline run successful"
```

---

## Notes

- Apptainer is designed for secure, rootless container use on HPC.
- `docker://` lets you run or build from Docker Hub images.
- Store `.sif` images in your `$HOME` or `$WORK` for persistent use.

---

## TL;DR Quick Start

```bash
module load apptainer/1.4.0
apptainer run docker://alpine echo "Hello from container"
```

---

## Resources

- Apptainer docs: https://apptainer.org
- VT ARC: https://arc.vt.edu
