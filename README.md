# Experiments-for-Tensor-Decomposition-with-Unaligned-Observations

In this repository, we provide the code scripts used to create the results in the following paper: 

Runshi Tang, Tamara Kolda, and Anru R. Zhang. Tensor Decomposition with Unaligned Observations

# Reproducing Experiments

## Setup

Requires R (only used to generate simulation data), Python3.11, and Python libraries imported in the scripts. 
The working directories and file paths in the scripts need to be adjusted accordingly to run the scripts. 
The experiments are performed on Winstat for big jobs and Slurm of the Social Science Computing Cooperative of UW-Madison. 

## Contents

### Simulation Study

* All the files are in the directory `simulation/`. 
* `generate5.R` is used to generate the data used for this simulation study. 
* Directory `rkhs` includes the files used for the simulation study involving RKHS-TD and S-RKHS-TD in the paper.
  - File `helper.py` contains the main function which implements RKHS-TD and S-RKHS-TD.
  - Directory `setting1` and `setting2` contain the codes for simulation studies of RKHS-TD and S-RKHS-TD with target rank 2 and 1 respectively, including performing the decomposition and plotting the figures in the paper. 
* Directory `grkhs_poisson` includes the files used for the simulation study involving GRKHS-TD and S-GRKHS-TD in the paper.
  - File `helper_counts_v3.py` contains the main function which implements GRKHS-TD and S-GRKHS-TD with Poisson loss.
  - File `run.py` performs the decomposition and `plot.py` plots the figures showing loss and time in the paper.

### Real Data Experiment

* All the files are in the directory `real_data/ECAM`.
* `microbiome.json` and `microbiome_relative_abundance.json` are ECAM data processed by CLR and relative abundance respectively.
* 










