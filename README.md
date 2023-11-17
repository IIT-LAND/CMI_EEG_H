# CMI_EEG_H
All code for Bertelsen et al., 2023, medRXiv
---
/n
## Directories to set up:
data: 
- pwd: ~/[root]/cmi_eeg_rs_H/data/tidy
/n
code 
- pwd: ~/[root]/cmi_eeg_rs_H/code
/n
results:
- pwd: ~/[root]/cmi_eeg_rs_H/results/reval/global
/n
plots: 
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsopen/adjH
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsopen/d
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsopen/pca
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsopen/pheno
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsopen/pls
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsopen/sigclust
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsopen/topoplots
/n
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsclosed/adjH
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsclosed/d
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsclosed/pheno
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsclosed/pca
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsclosed/pls
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsclosed/sigclust
- pwd: ~/[root]/cmi_eeg_rs_H/reval/global/rsclosed/topoplots
---
/n
## Pipeline for the analyses:
/n
### 01 run reval to identify subtypes
code:
- s01_reval_adjH.ipynb
/n
data to input:
- tidy_cov_adjusted_H.csv 
/n
result files (outputted to data/tidy directory):
- tidy_adjH_revalsubgroups_global.csv (main)
- tidy_adjH_revalsubgroups_rsopen_umap+cl.csv.csv (umap resting state eyes open)
- tidy_adjH_revalsubgroups_rsclosed_umap+cl.csv.csv (umap resting state eyes closed)
/n
* note: if you have a single file with H+pheno data you'll have to adjust this script
---

### 02 check cluster consistency across rs eye conditions
code:
- s02_reval_adjH_subtypes_consistency.Rmd

data to input:
- tidy_adjH_revalsubgroups_global.csv

plots:
- plots/reval/global/02_subtype_consistency_alluvialplot.pdf
---

### 03 check subtypes with SigClust
code:
- s03_reval_adjH_sigClust_subtypes_validation.Rmd

data to input:
- tidy_H_revalsubgroups_rsopne_umap+cl.csv
- tidy_H_revalsubgroups_rsclosed_umap+cl.csv

result files:
- 03_sigclust10000x_results.csv

plots (same for rsclosed):
- plots/reval/global/rsopen/sigclust/03_adjH_UMAP_training.pdf
- plots/reval/global/rsopen/sigclust/03_adjH_UMAP_validation.pdf
---

### 04 run linear model at the electrode level 
code:
- cmi_02_H_reval_G30nC2_asd_males_channel_subtypeAnalysis_31.01.2023.Rmd

data to input:
- tidy_adjH_revalsubgroups_global.csv

result files:
- rsopen_electrodes_subtype*age_06.Mar.2023.csv
- rsclosed_electrodes_subtype*age_06.Mar.2023.csv


plots (same for rsclosed):
H values plots:
- plots/reval/global/rsopen/adjH/04_adjH_byblock_%s.pdf
- plots/reval/global/rsopen/adjH/04_adjH_blockscollapsed_%s.pdf
- plots/reval/global/rsopen/adjH/04_adjH_age_%s.pdf
- plots/reval/global/rsopen/adjH/04_adjH_age2_%s.pdf
- plots/reval/global/rsopen/adjH/04_adjH_age_allElectrodes_facetBlocks.png
---

### 05 plotting topoplots: average H, mean abs H differences, F and p-val
code:
- s05_reval_adjH_electrode_subtypes_topoplots.mlx

data to input:
- tidy_adjH_revalsubgroups_global.csv
- 04_adjH_rsopen_electrodes_subtype*age.csv
- 04_adjH_rsclosed_electrodes_subtype*age.csv

plots (same for rsclosed):
- plots/reval/global/rsopen/topoplots/05_adjH_topoplot_mean_H_td.jpg
- plots/reval/global/rsopen/topoplots/05_adjH_topoplot_mean_H_a1.jpg
- plots/reval/global/rsopen/topoplots/05_adjH_topoplot_mean_H_a2.jpg
- plots/reval/global/rsclosed/topoplots/05_adjH_topoplot_mean_H_abs_diff_td-a1.jpg
- plots/reval/global/rsclosed/topoplots/05_adjH_topoplot_mean_H_abs_diff_td-a2.jpg
- plots/reval/global/rsclosed/topoplots/05_adjH_topoplot_lme_age_F.jpg
- plots/reval/global/rsclosed/topoplots/05_adjH_topoplot_lme_age_pval_fdr.jpg
---

### 06 run PCA 
code:
- s06_reval_adjH_runPCA.mlx

data to input:
- tidy_H_reval_asd_males_revalsubgroups_global_01.02.2023.csv


result files:
- data/tidy/tidy_H_reval_asd_males_PCA_rsopen_02.Feb.2023.csv
- data/tidytidy_H_reval_asd_males_PCA_rsclosed_02.Feb.2023.csv
- results/reval/global/rsopen_PCA_variance_explained.xlsx
- results/reval/global/rsclosed_PCA_variance_explained.xlsx

plots (scree plots, loadings over pc, scores by electrode topoplots):
- plots/reval/global/rsopen/topoplots/rsopen_H_PCA_pct_variance_explained.jpg
- plots/reval/global/rsopen/topoplots/rsopen_H_PCA_loadings_over_pc.jpg
- plots/reval/global/rsopen/topoplots/rsopen_H_PCA_topoplots.jpg
- plots/reval/global/rsopen/topoplots/topoplot_H_PC1.jpg
- etc

---

### plot effect size differences by subtype and across the five blocks for electrodes (all 93) and PCs (1, 2 and 3)

code:
- cmi_04_H_reval_G30nC2_asd_males_effectsize_comparison_plots.Rmd
- cmi_04_H_reval_G30nC2_asd_males_effectsize_comparison_plots.html


data to input:
- rsopen_electrodes_subtype*age_06.Mar.2023.csv
- rsclosed_electrodes_subtype*age_06.Mar.2023.csv
- rsopen_PCA_subtype*age_24.Apr.2023.csv
- rsclosed_PCA_subtype*age_24.Apr.2023.csv


plots:
- plots/reval/global/rsopen/d/H_rsopen_d_byelectrode_acrossblocks_2col.pdf
- plots/reval/global/rsclosed/d/H_rsclosed_d_byelectrode_acrossblocks_2col.pdf
- plots/reval/global/rsopen/d/pca_pc1_rsopen_d_acrossblocks_2col.pdf
- plots/reval/global/rsclosed/d/pca_pc1_rsclosed_d_acrossblocks_2col.pdf.pdf

---

### plot electrode effects sizes in the topoplots

code:
- cmi_04_H_reval_G30nC2_asd_males_topoplots_effectsizes_16032023.mlx

data to input:
- rsopen_electrodes_subtype*age_06.Mar.2023.csv
- rsclosed_electrodes_subtype*age_06.Mar.2023.csv

plots:
- plots/reval/global/rsopen/topoplots/topoplot_tdva1_mean_d_2col.jpg
- plots/reval/global/rsopen/topoplots/topoplot_tdva2_mean_d_2col.jpg
- plots/reval/global/rsclosed/topoplots/topoplot_tdva1_mean_d_2col.jpg
- plots/reval/global/rsclosed/topoplots/topoplot_tdva2_mean_d_2col.jpg


---


### run linear model at the PCA level

code:
- cmi_05_H_reval_asd_males_PCA_subtypeAnalyses_31.01.2023.Rmd
- cmi_05_H_reval_asd_males_PCA_subtypeAnalyses_31.01.2023.html

data to input:
- tidy_H_reval_asd_males_PCA_rsopen_02.Feb.2023.csv
- tidy_H_reval_asd_males_PCA_rsclosed_02.Feb.2023.csv

result files:
- rsopen_PCA_subtype*age_24.Apr.2023.csv
- rsclosed_PCA_subtype*age_24.Apr.2023.csv

---

### check for pheno differences between the subtypes

code:
- cmi_06_H_reval_asd_males_phenoAnalyses_31.01.2023.Rmd
- cmi_06_H_reval_asd_males_phenoAnalyses_31.01.2023.Rmd

data to input:
- tidy_H_reval_asd_males_revalsubgroups_global_01.02.2023.csv

result files:
- restingstate_pheno_summary.csv
- restingstate_pheno_analyses.csv

some plots (age, fiq, srs, rbs, scq, assq):
- plots/reval/global/rsopen/pheno/age.pdf

  
---

### extended pheno plotting

code:
- cmi_06_H_reval_G30nC2_asd_males_phenoAnalyses_Extended_25.02.2023.Rmd
- cmi_06_H_reval_G30nC2_asd_males_phenoAnalyses_Extended_25.02.2023.html

data to input:
- tidy_H_reval_asd_males_revalsubgroups_global_01.02.2023.csv

plots:
- plots/reval/global/rsopen/pheno_extended/c3sr_ab_T.jpg

---
