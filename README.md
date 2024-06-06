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

## Directories to set up:
data: 
- pwd: ~/[root]/cmi_eeg_rs_H/data/tidy  

code 
- pwd: ~/[root]/cmi_eeg_rs_H/code
 
results:
- pwd: ~/[root]/cmi_eeg_rs_H/results/reval/global  

plots: 
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsopen/adjH
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsopen/d
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsopen/pca
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsopen/pheno
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsopen/pls
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsopen/sigclust
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsopen/topoplots
  - ...
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsclosed/adjH
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsclosed/d
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsclosed/pheno
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsclosed/pca
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsclosed/pls
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsclosed/sigclust
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsclosed/topoplots

---  

## In-silico computational modeling

`_insilico_analysis.Rmd` or `_insilico_analysis.html` 

## In-vivo chemogenetic analyses

`_invivo_analysis.Rmd` or `_invivo_analysis.html` 


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

data to input:
- `tidy_cov_adjusted_H.csv` 

result files (outputted to `data/tidy directory`):
- `tidy_adjH_revalsubgroups_global.csv` (main)
- `tidy_adjH_revalsubgroups_rsopen_umap+cl.csv` (umap resting state eyes open)
- `tidy_adjH_revalsubgroups_rsclosed_umap+cl.csv` (umap resting state eyes closed)

* note: if you have a single file with H+pheno data you'll have to adjust this script

---

### 02 check cluster consistency across rs eye conditions
code:
- `s02_reval_adjH_subtypes_consistency.Rmd`

data to input:
- `tidy_adjH_revalsubgroups_global.csv`

plots:
- `plots/reval/global/02_subtype_consistency_alluvialplot.pdf`

---

### 03 check subtypes with SigClust
code:
- `s03_reval_adjH_sigClust_subtypes_validation.Rmd`

data to input:
- `tidy_H_revalsubgroups_rsopen_umap+cl.csv`
- `tidy_H_revalsubgroups_rsclosed_umap+cl.csv`

result files:
- `03_sigclust10000x_results.csv`

plots (same for rsclosed):
- `plots/reval/global/rsopen/sigclust/03_adjH_UMAP_training.pdf`
- `plots/reval/global/rsopen/sigclust/03_adjH_UMAP_validation.pdf`

---

### 04 run linear model at the electrode level 
code:
- `s04_reval_adjH_electrode_subtypes_analysis.Rmd`

data to input:
- `tidy_adjH_revalsubgroups_global.csv`

result files:
- `rsopen_electrodes_subtype*age_06.Mar.2023.csv`
- `rsclosed_electrodes_subtype*age_06.Mar.2023.csv`


plots (same for rsclosed):
- `plots/reval/global/rsopen/adjH/04_adjH_byblock_%s.pdf`
- `plots/reval/global/rsopen/adjH/04_adjH_blockscollapsed_%s.pdf`
- `plots/reval/global/rsopen/adjH/04_adjH_age_%s.pdf`
- `plots/reval/global/rsopen/adjH/04_adjH_age2_%s.pdf`
- `plots/reval/global/rsopen/adjH/04_adjH_age_allElectrodes_facetBlocks.png`

---

### 05 plotting topoplots: average H, mean abs H differences, F and p-val
code:
- `s05_reval_adjH_electrode_subtypes_topoplots.mlx`

data to input:
- `tidy_adjH_revalsubgroups_global.csv`
- `04_adjH_rsopen_electrodes_subtype*age.csv`
- `04_adjH_rsclosed_electrodes_subtype*age.csv`

plots (same for rsclosed):
- `plots/reval/global/rsopen/topoplots/05_adjH_topoplot_mean_H_td.jpg`
- `plots/reval/global/rsopen/topoplots/05_adjH_topoplot_mean_H_a1.jpg`
- `plots/reval/global/rsopen/topoplots/05_adjH_topoplot_mean_H_a2.jpg`
- `plots/reval/global/rsclosed/topoplots/05_adjH_topoplot_mean_H_abs_diff_td-a1.jpg`
- `plots/reval/global/rsclosed/topoplots/05_adjH_topoplot_mean_H_abs_diff_td-a2.jpg`
- `plots/reval/global/rsclosed/topoplots/05_adjH_topoplot_lme_age_F.jpg`
- `plots/reval/global/rsclosed/topoplots/05_adjH_topoplot_lme_age_pval_fdr.jpg`

---

### 06 run PCA and reconstruct H data from PC4
code:
- `s06_reval_adjH_runPCA.mlx`

data to input:
- `tidy_adjH_revalsubgroups_global.csv`

result files:
- `data/tidy/tidy_adjH_rsopen_PCA.csv`
- `data/tidy/tidy_adjH_rsclosed_PCA.csv`
- `data/tidy/tidy_adjH_rsopen_reconstructedHfromPC4.csv`

plots (same for rsclosed):
- `plots/reval/global/rsopen/pca/06_adjH_PCA_pct_variance_explained.jpg`
- `plots/reval/global/rsopen/pca/06_adjH_PCA_cum_pct_variance_explained.jpg`
- `plots/reval/global/rsopen/pca/06_adjH_PCA_loadings_over_PCs.jpg`
- `plots/reval/global/rsopen/topoplots/06_adjH_topoplots_PCs1-8.jpg`
- `plots/reval/global/rsopen/topoplots/06_adjH_topoplots_PC1.jpg`
- `plots/reval/global/rsopen/topoplots/06_adjH_topoplots_PC2.jpg`
- `plots/reval/global/rsopen/topoplots/06_adjH_topoplots_PC3.jpg`
- `plots/reval/global/rsopen/topoplots/06_adjH_topoplots_PC4.jpg`

---

### 07 run linear model at the PCA level
code:
- `s07_reval_adjH_PCA_subtypes_analysis.Rmd`

data to input:
- `tidy_adjH_rsopen_PCA.csv`
- `tidy_adjH_rsclosed_PCA.csv`

result files:
- `07_adjH_rsopen_PCA_subtype*age.csv`
- `07_adjH_rsclosed_PCA_subtype*age.csv`

plots (same for rsclosed):
- `plots/reval/global/rsopen/pca/07_adjH_PCAscores_%s.pdf`
- `plots/reval/global/rsopen/pca/07_adjH_PCAscores_%s_blockscollapsed.pdf`
- `plots/reval/global/rsopen/pca/07_adjH_PCAscores_%s_age_groupscollapsed.pdf`
- `plots/reval/global/rsopen/pca/07_adjH_PCAscores_%s_age_blockscollapsed.pdf`

---

### 08 plot effect size differences by subtype and across the five blocks for electrodes (all 93) and PCs (1, 2 and 3)
code:
- `s08_reval_adjH_effectsizes_plots.Rmd`

data to input:
- `rsopen_electrodes_subtype*age_06.Mar.2023.csv`
- `rsclosed_electrodes_subtype*age_06.Mar.2023.csv`
- `rsopen_PCA_subtype*age_24.Apr.2023.csv`
- `rsclosed_PCA_subtype*age_24.Apr.2023.csv`

plots:
- `plots/reval/global/rsopen/d/H_rsopen_d_byelectrode_acrossblocks_2col.pdf`
- `plots/reval/global/rsclosed/d/H_rsclosed_d_byelectrode_acrossblocks_2col.pdf`
- `plots/reval/global/rsopen/d/pca_pc1_rsopen_d_acrossblocks_2col.pdf`
- `plots/reval/global/rsclosed/d/pca_pc1_rsclosed_d_acrossblocks_2col.pdf`

---

### 09 calculate and plot effect sizes for rec H data (PC1+PC4)
code: 
- `s09_reval_adjH_effectsizes_recH.Rmd`

data to input: 
- `tidy_adjH_rsopen_reconstructedHfromPC1&4.csv`

results:
- `09_adjH_rsopen_Hrec-PC1PC4_effectsizes.csv`

plots:
- `plots/reval/global/rsopen/d/09_adjH_d_Hrec-PC1PC4_acrossblocks.pdf`

---

### 10 topoplots for effects sizes (H + recH data)
code: 
- `s10_reval_adjH_effectsizes_topoplots.mlx`

data to input:
- `04_adjH_rsopen_electrodes_subtype*age.csv`
- `04_adjH_rsopen_electrodes_subtype*age.csv`
- `09_adjH_rsopen_Hrec-PC1PC4_effectsizes.csv`

results:

plots (same for rsclosed):
- `plots/reval/global/rsopen/d/10_adjH_topoplot_mean_d_electH_tdva1.jpg`
- `plots/reval/global/rsopen/d/10_adjH_topoplot_mean_d_electH_tdva2.jpg`
- `plots/reval/global/rsopen/d/10_adjH_topoplot_mean_d_PC1&4recH_tdva1.jpg`
- `plots/reval/global/rsopen/d/10_adjH_topoplot_mean_d_PC1&4recH_tdva2.jpg`

---

### 11 Pheno analyses
code: 
- `s11_reval_adjH_pheno_analysis.Rmd`

data to input: 
- `tidy_adjH_revalsubgroups_global.csv`

results:
- `11_adjH_rsopen_pheno_analyses.csv`
- `11_adjH_rsclosed_pheno_analyses.csv`
- `11_adjH_restingstate_pheno_summary.csv`

plots (same for rsclosed):
- `plots/reval/global/rsopen/pheno/11_adjH_age.pdf`
- `plots/reval/global/rsopen/pheno/11_adjH_fiq.pdf`
- `plots/reval/global/rsopen/pheno/11_adjH_assq_total.pdf`
- `plots/reval/global/rsopen/pheno/11_adjH_rbs_total.pdf`
- `plots/reval/global/rsopen/pheno/11_adjH_scq_total.pdf`
- `plots/reval/global/rsopen/pheno/11_adjH_srs_total_T.pdf`

---

### 12 prepare dataframes for PLS
code: 
- `s12_reval_adjH_prepdata4PLS.Rmd`

data to input:
- `tidy_adjH_revalsubgroups_global.csv`
- `mvl_tidy_pheno_4pls.csv`

results (same for rsclosed):
- `tidy_rsopen_pheno_4pls.csv`
- `wide_rsopen_adjH_b1.csv`
- `wide_rsopen_adjH_b1_td.csv`
- `wide_rsopen_adjH_b2_td.csv`
- `wide_rsopen_adjH_b3_td.csv`
- `wide_rsopen_adjH_b4_td.csv`
- `wide_rsopen_adjH_b5_td.csv`
- `wide_rsopen_adjH_b1_a1.csv`
- `wide_rsopen_adjH_b2_a1.csv`
- `wide_rsopen_adjH_b3_a1.csv`
- `wide_rsopen_adjH_b4_a1.csv`
- `wide_rsopen_adjH_b5_a1.csv`
- `wide_rsopen_adjH_b1_a2.csv`
- `wide_rsopen_adjH_b2_a2.csv`
- `wide_rsopen_adjH_b3_a2.csv`
- `wide_rsopen_adjH_b4_a2.csv`
- `wide_rsopen_adjH_b5_a2.csv`

---

### 13, 14 run PLS
code: 
- `s13_reval_adjH_runPLS.m`
- `s14_reval_adjH_compute_bootstrap_ci_pls.m`

data to input (same for rsclosed):
- `tidy_rsopen_pheno_4pls.csv`
- `wide_rsopen_adjH_b1.csv`
- `wide_rsopen_adjH_b1_td.csv`
- `wide_rsopen_adjH_b2_td.csv`
- `wide_rsopen_adjH_b3_td.csv`
- `wide_rsopen_adjH_b4_td.csv`
- `wide_rsopen_adjH_b5_td.csv`
- `wide_rsopen_adjH_b1_a1.csv`
- `wide_rsopen_adjH_b2_a1.csv`
- `wide_rsopen_adjH_b3_a1.csv`
- `wide_rsopen_adjH_b4_a1.csv`
- `wide_rsopen_adjH_b5_a1.csv`
- `wide_rsopen_adjH_b1_a2.csv`
- `wide_rsopen_adjH_b2_a2.csv`
- `wide_rsopen_adjH_b3_a2.csv`
- `wide_rsopen_adjH_b4_a2.csv`
- `wide_rsopen_adjH_b5_a2.csv`

results:
- `13_pls_rsopen_ALL_BSR_LV1.csv`
- `13_pls_rsopen_ALL_BSRrev_LV1.csv`
- `13_pls_rsopen_ALL.mat`
- `13_pls_rsclosed_ALL_BSR_LV1.csv`
- `13_pls_rsclosed_ALL_BSRrev_LV1.csv`
- `13_pls_rsclosed_ALL.mat`
- `14_rsopen_ALL_bootlim_data4plotting_LV1_ci95.csv`
- `14_rsopen_ALLrev_bootlim_data4plotting_LV1_ci95.csv`
- `14_rsclosed_ALL_bootlim_data4plotting_LV1_ci95.csv`
- `14_rsclosed_ALLrev_bootlim_data4plotting_LV1_ci95.csv`

---

### 16
code: 
- `s16_reval_adjH_plot_pls_results.Rmd`

data to input:
- `13_pls_rsopen_ALL_BSRrev_LV1.csv`
- `13_pls_rsclosed_ALL_BSR_LV1.csv`
- `14_rsopen_ALLrev_bootlim_data4plotting_LV1_ci95.csv`

plots (same for rsclosed):
- `plots/reval/global/rsopen/pls/16_rsopen_ALL_BSR_LV1_corr.pdf`
- `plots/reval/global/rsopen/pls/16_rsclosed_ALL_BSR_LV1_corr.pdf`
- `plots/reval/global/rsopen/pls/16_rsopen_LV1_corr.pdf`
- `plots/reval/global/rsopen/pls/16_rsclosed_LV1_corr.pdf`
- `plots/reval/global/rsopen/pls/16_rsopen_LV1_corr.pdf`
- `plots/reval/global/rsopen/pls/16_rsclosed_LV1_corr.pdf`

---

### 17
code: 
- `s17_reval_adjH_pls_topoplots.mlx`

data to input:
- `13_pls_rsopen_ALL_BSRrev_LV1.csv`
- `13_pls_rsclosed_ALL_BSR_LV1.csv`

plots (same for rsclosed):
- `plots/reval/global/rsopen/pls/17_adjH_topoplot_PLS_BSR.jpg`
- `plots/reval/global/rsopen/pls/17_adjH_corr_matrix_byblock_BSR.jpg`

---
