#!/bin/bash

rootdir=/media/DATA/RAW/cmihbn
codedir=${rootdir}/code

# Step 1: Download and Cleaning ------------------------------------------------
cd ${codedir}

# Download
python _01_cmi_download.py server TD batch 1 10
python _01_cmi_download.py server ASD batch 1 10

# Cleaning
python _01_cmi_cleaning.py server step1 batch
python _01_cmi_cleaning.py server step2 batch
python _01_cmi_cleaning.py server step3 batch

# Raw plots
python _01_cmi_rawplots.py server batch

# Step 2: Preprocessing --------------------------------------------------------

# Preprocessing
python _02_cmi_preprocessing.py server batch step1

# *** Manual QC step - look through preproc reports and make tidy_manual_qc.csv ***

# Step 3: Postprocessing --------------------------------------------------------
# Postprocessing
python _03_cmi_postprocessing_hurst.py server batch

# Step 4: Tidying data ---------------------------------------------------------
python _04_cmi_tidy_dataquality.py server batch
python _04_cmi_tidy_hurst.py server batch
