# Electrophysiologically-defined excitation-inhibition autism neurosubtypes 

This repository has all the code and tidy data for the analyses in **Bertelsen, N., Mancini, G., Sastre-Yagüe, D., Vitale, A., Lorenz, G. M., Malerba, S. B., Bolis, D., Mandelli, V., Martínez-Cañada, P., Gozzi, A., Panzeri, S., & Lombardo, M. V. Electrophysiologically-defined excitation-inhibition autism neurosubtypes. medRxiv. doi:10.1101/2023.11.22.23298729.** (https://www.medrxiv.org/content/10.1101/2023.11.22.23298729v1)

The data utilized in this work comes from the publicly available CMI-HBN dataset (http://fcon_1000.projects.nitrc.org/indi/cmi_healthy_brain_network/). For more details about the CMI-HBN data see Alexander et al., 2017, Scientific Data (https://www.nature.com/articles/sdata2017181).

The code directory within this repository is organized to have all the main initial steps (within the `code/pp` directory) and will also have all the other code for downstream data analysis within it. Below is a description of each step of the analyses, with the code, data, and results to expect/use.

---

## General requirements, assumptions, dependencies

Analyses heavily depend on use of `MATLAB`, `R`, and `Python`. For example, `reval` is implemented in `Python`. `MATLAB` is heavily relied upon for EEG preprocessing, computation of the Hurst exponent (H), and PLS analysis. `R` is heavily relied upon for main downstream statistical analysis (e.g., linear mixed effect modeling, plotting data, SigClust analysis). Below is a list of primary libraries/toolboxes that primary components of the data analysis heavily rely upon.  

  + **EEGLAB** (https://eeglab.org) - Used ubiquitously throughout for processing and handling EEG data.

  + **reval** (https://github.com/IIT-LAND/reval_clustering) - Used for stability-based relative clustering validation analyses.

  + **nonfractal** (https://github.com/wonsang/nonfractal) - Used for computation of H.

  + **PLS toolbox in MATLAB** (https://www.rotman-baycrest.on.ca/index.php?section=84) - Used for the PLS analysis.

---

## In-silico computational modeling

`_insilico_analysis.Rmd` or `_insilico_analysis.html` 

---

## In-vivo chemogenetic analyses

`_invivo_analysis.Rmd` or `_invivo_analysis.html` 

---

## Pipeline for the human analyses:

### **pp** download, clean, preprocess, and postproc H
code:
- code for running these steps is within the `code/pp` directory

`_00_master.sh` specifies sequence of steps for running code to do all steps
`_0*.py` scripts implement each step 

---

### 01 run reval to identify subtypes
code:
- `s01_reval_adjH.ipynb`

---

### 02 check cluster consistency across rs eye conditions

code:
- `s02_reval_adjH_subtypes_consistency.Rmd`

---

### 03 check subtypes with SigClust

code:
- `s03_reval_adjH_sigClust_subtypes_validation.Rmd`

---

### 04 run linear model at the electrode level 

code:
- `s04_reval_adjH_electrode_subtypes_analysis.Rmd`

---

### 05 plotting topoplots: average H, mean abs H differences, F and p-val

code:
- `s05_reval_adjH_electrode_subtypes_topoplots.mlx`

---

### 06 run PCA and reconstruct H data from PC4

code:
- `s06_reval_adjH_runPCA.mlx`

---

### 07 run linear model at the PCA level

code:
- `s07_reval_adjH_PCA_subtypes_analysis.Rmd`

---

### 08 plot effect size differences by subtype and across the five blocks for electrodes (all 93) and PCs (1, 2 and 3)

code:
- `s08_reval_adjH_effectsizes_plots.Rmd`

---

### 09 calculate and plot effect sizes for rec H data (PC1+PC4)

code: 
- `s09_reval_adjH_effectsizes_recH.Rmd`

---

### 10 topoplots for effects sizes (H + recH data)

code: 
- `s10_reval_adjH_effectsizes_topoplots.mlx`

---

### 11 Pheno analyses

code: 
- `s11_reval_adjH_pheno_analysis.Rmd`

---

### 12 prepare dataframes for PLS

code: 
- `s12_reval_adjH_prepdata4PLS.Rmd`

---

### 13, 14 run PLS

code: 
- `s13_reval_adjH_runPLS.m`
- `s14_reval_adjH_compute_bootstrap_ci_pls.m`

---

### 16

code: 
- `s16_reval_adjH_plot_pls_results.Rmd`

---

### 17

code: 
- `s17_reval_adjH_pls_topoplots.mlx`

---
