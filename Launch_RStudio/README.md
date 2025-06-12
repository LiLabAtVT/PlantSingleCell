
# Launching RStudio via ARC Open OnDemand (Virginia Tech)

This guide walks you through the steps required to access RStudio on the ARC Open OnDemand portal at Virginia Tech.

---

## Step 1: Connect to VPN

Use the **Ivanti Secure Access Client** to connect to the VT VPN.

<!-- Original Markdown -->
<!-- ![Step 1 - VPN Connection](VPN.png) -->

<!-- Updated with smaller size -->
<img src="VPN.png" alt="Step 1 - VPN Connection" width="400"/>

---

## Step 2: Login with VT Credentials

When prompted, enter your **VT PID** as the username and your **VT password**.

<!-- Original Markdown -->
<!-- ![Step 2 - VT Credentials](Password.png) -->

<!-- Updated with smaller size -->
<img src="Password.png" alt="Step 2 - VT Credentials" width="400"/>

---

## Step 3: Open ARC Open OnDemand

Navigate to the following URL in your browser:

[https://ood.arc.vt.edu](https://ood.arc.vt.edu)

This will bring you to the ARC Open OnDemand dashboard.

<!-- Original Markdown -->
<!-- ![Step 3 - ARC Dashboard](RStudio_app.png) -->

<!-- Updated with smaller size -->
<img src="RStudio_app.png" alt="Step 3 - ARC Dashboard" width="400"/>

---

## Step 4: Launch RStudio

From the Open OnDemand dashboard, locate the **RStudio Server** application.

![Step 4 - Launch RStudio](RStudio_launch.png)

### 4.1: Select Cluster

You can choose from three available clusters:

- **TinkerCliffs**: CPU-based (general use)
- **Owl**: CPU-based 
- **Falcon**: GPU-based 



### 4.2: Partition

Partitions define the scheduling queues available for a given cluster. Choose the appropriate one based on your selected cluster
Refer to ARC documentation for details:  
[VT ARC Compute Resources](https://www.docs.arc.vt.edu/resources/compute.html)



### 4.3: Account

This refers to your **allocation name**. For this guide, we use:

```text
plantsinglecell
```

You can view your active allocations with:

```bash
quota
```



### 4.4: Number of Hours

This sets how long your RStudio session will last. You can choose between **1 to 48 hours**.  


### 4.5: R Version

Choose the appropriate R version available on the cluster.


### 4.6: Number of Cores per Node

This defines the number of **CPU cores** allocated for your session.  

