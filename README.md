# Experiments for Tensor Decomposition with Unaligned Observations

This repository contains the code scripts used in the research presented in the paper:

- **Title:** Tensor Decomposition with Unaligned Observations
- **Authors:** Runshi Tang, Tamara Kolda, and Anru R. Zhang

# Reproducing Experiments

## Setup
The experiments require:
- R (for generating simulation data and implementing CP Decomposition and Tabular FTSVD)
- Python 3.11
- Python libraries: `os`, `numpy`, `json`, `time`, `scipy`, `random`, `pandas`, `sys`, `glob`, `matplotlib`

Users need to adjust working directories and file paths in the scripts to match their environments. The experiments were originally conducted in Winstat for Big Jobs and the Slurm of the Social Science Computing Cooperative at the University of Wisconsin-Madison.

## Contents

* `algorithms.py` contains all algorithms discussed in the paper:
    - `RKHS` implements RKHS-TD (with `sketch = False`) and S-RKHS-TD with `sketch = True`). The parameter `lambda_parameter = x` sets the penalty parameter $\lambda = x$ in the paper.  
    - `RKHS_generalized_counts` implements GRKHS-TD with Poisson loss.
    - `RKHS_generalized_sketch_counts` implements S-GRKHS-TD with Poisson loss.
    - `RKHS_generalized_sketch_beta_divergnece` implements S-GRKHS-TD with beta divergnece loss.
    - `RKHS_generalized_beta_divergnece` implements GRKHS-TD with beta divergence loss.
    - In `RKHS_generalized_...`, the parameter `lambda_parameter = x` sets the parameter $C_\xi = x$ in the feasible set.

### Simulation Study

- The `simulation/` directory contains all related files.
- `rkhs/`: Includes files for the simulations on RKHS-TD and S-RKHS-TD.
  - `generate7.R`: Generates simulation data.
  - File `helpers.py` Implements RKHS-TD and S-RKHS-TD. Setting `sketch = False` implements RKHS-TD; setting `sketch = True` implements S-RKHS-TD. 
  - `setting1/`: Contains codes for simulation with RKHS-TD, including decomposition (`run.py`) and plotting results (`plot.py`).
- `grkhs_poisson/`: Contains files for GRKHS-TD and S-GRKHS-TD simulations.
  - `generate8.R`: Generates simulation data.
  - `helper_counts_v2.py`: Implements GRKHS-TD and S-GRKHS-TD with Poisson loss.
  - `run.py`: Performs decomposition, and `plot.py` plots loss and time figures.
- To replicate the simulation results in the paper, first run `generate7.R` or `generate8.R` to generate the simulated data. Next, run `run.py` to generate and save the simulation results. Finally, run `plot.py` for figures. 
  - Note that the first few lines in each code script need to be changed accordingly to set working directory to the location of the code script. 

### Real Data Experiment

- Located in `real_data/ECAM/`.
- `microbiome.json` and `microbiome_relative_abundance.json`: Processed ECAM data by centered-log-ratio and relative abundance, respectively.
- Directory `rkhs/`: Contains files for the real data study with RKHS-TD and S-RKHS-TD.
  - File `helpers.py` contains the main function which implements RKHS-TD and S-RKHS-TD. Setting `sketch = False` will implement RKHS-TD, and `sketch = True` will implement S-RKHS-TD. 
  - Directory `rx/s2_xx` contains the codes for RKHS-TD and S-RKHS-TD with target rank r = x and s2 = xx respectively, including performing the decomposition (`run.py`) and plotting the figures in the paper (`plot.py`). 
- Directory `grkhs_beta_d/`: Includes files for the real data study with GRKHS-TD and S-GRKHS-TD with Beta divergence loss.
  - `helper_betad.py`: Implements GRKHS-TD and S-GRKHS-TD with Beta divergence loss.
  - `rx/s2_xx`: Contains codes for GRKHS-TD and S-GRKHS-TD with target rank `r = x` and `s2 = xx`, including decomposition (`run.py`) and plotting (`plot.py`).
- The above two directories contain the real data experiments of Figure 5 and Figure 6 in the paper. To replicate, run each `run.py` to compute and save the results, and then run `plot.py` for figures. 
  - Note that the first few lines in each code script need to be changed accordingly to set working directory to the location of the code script. 
- Directory `compare/`: Files for comparing CP Decomposition, Tabular FTSVD, GCP, SRKHS−TD, and G-SRKHS−TD.
  - `helpers_v2.py`: Implements S-RKHS-TD.
  - `helper_betad.py`: Implements S-GRKHS-TD with Beta divergence loss.
  - `RKHS.py`: Performs SRKHS−TD on the ECAM data. 
  - `GRKHS.py`: Performs G-SRKHS−TD on the ECAM data. 
  - `FTSVD.R`: Contains the main function which implements Tabular FTSVD.
  - `realdata2.R`: Implements Tabular FTSVD and CP Decomposition. This file and `FTSVD.R` are sourced from [`https://github.com/Rungang/functional_tensor_svd`](https://github.com/Rungang/functional_tensor_svd) and [`https://cran.r-project.org/web/packages/rTensor/`](https://cran.r-project.org/web/packages/rTensor/).
  - `GCP.m`: Implements GCP with Beta divergence loss. The code is sourced from [`https://www.tensortoolbox.org/gcp_opt_doc.html`](https://www.tensortoolbox.org/gcp_opt_doc.html). 
  - Some `.csv` and `.json`: Contains raw ECAM data (used by `realdata2.R`) and more subject information. 
  - `2.R`: Generates results in the paper (require results from `realdata2.R`, `GCP.m`, `RKHS.py` and `GRKHS.py`).
  - This directory contains the real data experiments of Figure 7 and Figure 8 in the paper. To replicate, run `realdata2.R`, `GCP.m`, `RKHS.py` and `GRKHS.py`, and then run `2.R` for Figures. 
  - Note that the first few lines in each code script need to be changed accordingly to set working directory to the location of the code script. 
  - We keep all the intermediate results in this repo, so the user can directly run `2.R` to generate the Figures in the paper. 
