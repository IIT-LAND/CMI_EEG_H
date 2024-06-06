#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# grab start and end data release input arguments
group2use = as.character(args[1])
start_data_release = as.numeric(args[2])
end_data_release = as.numeric(args[3])

datareleases2use = seq(start_data_release,end_data_release)
if (length(datareleases2use)==1){
  fstem = sprintf("R%d", start_data_release)
} else {
  fstem = sprintf("R%d_to_R%d", start_data_release, end_data_release)
} # if (length(datareleases2use)==1){


# load packages
library(easypackages)
libraries("here","tidyverse")
source(here("code","cmi_utils.R"))

# directories
code_dir = here("code")
pheno_dir = here("data","pheno")
tmp_dir = here("tmp")

# read in main phenotypic data file
pheno_file = "data-2022-07-01T06_27_26.878Z.csv"
main_pheno_file = file.path(pheno_dir,pheno_file)
main_pheno_data = read.csv(main_pheno_file)

# find subject IDs for specific diagnostic groups
if (group2use=="TD"){
  td_subs = dx_finder(pheno_data = main_pheno_data, dx2find = "No Diagnosis Given")
} else if (group2use=="ADHD"){
  adhd_subs = dx_finder(pheno_data = main_pheno_data, dx2find = "ADHD")
} else if (group2use=="ASD"){
  asd_subs = dx_finder(pheno_data = main_pheno_data, dx2find = "Autism")
} else if (group2use=="ID"){
  id_subs = dx_finder(pheno_data = main_pheno_data, dx2find = "Intellectual")
} else if (group2use=="LANG"){
  lang_subs = dx_finder(pheno_data = main_pheno_data, dx2find = "Language")
}

# find subjects that have EEG data
for (ifile in 1:length(datareleases2use)){

  data_release_num = datareleases2use[ifile]

  if (data_release_num<3){
    fname = file.path(pheno_dir, sprintf("HBN_R%d_1_Pheno.csv",data_release_num))
  } else {
    fname = file.path(pheno_dir, sprintf("HBN_R%d_Pheno.csv",data_release_num))
  } # if (data_release_num<3){

  tmp_data = read.csv(fname)
  if (ifile==1){
    eeg_pheno_data = tmp_data
  } else{
    eeg_pheno_data = rbind(eeg_pheno_data,tmp_data)
  } # if (ifile==1){

} # for (data_release_num in datareleases2use){

# write subject lists to file
colnames(eeg_pheno_data)[1] = "subid"
if (group2use=="TD"){
  td_eeg_pheno_data = eeg_pheno_data %>% filter(is.element(subid, td_subs$subid)) %>% distinct()
  write(as.character(td_eeg_pheno_data$subid), file = file.path(tmp_dir,sprintf("td_sublist_%s.csv",fstem)), sep = ",")
} else if (group2use=="ADHD"){
  adhd_eeg_pheno_data = eeg_pheno_data %>% filter(is.element(subid, adhd_subs$subid)) %>% distinct()
  write(as.character(adhd_eeg_pheno_data$subid), file = file.path(tmp_dir,sprintf("adhd_sublist_%s.csv",fstem)), sep = ",")
} else if (group2use=="ASD"){
  asd_eeg_pheno_data = eeg_pheno_data %>% filter(is.element(subid, asd_subs$subid)) %>% distinct()
  write(as.character(asd_eeg_pheno_data$subid), file = file.path(tmp_dir,sprintf("asd_sublist_%s.csv",fstem)), sep = ",")
} else if (group2use=="ID"){
  id_eeg_pheno_data = eeg_pheno_data %>% filter(is.element(subid, id_subs$subid)) %>% distinct()
  write(as.character(id_eeg_pheno_data$subid), file = file.path(tmp_dir,sprintf("id_sublist_%s.csv",fstem)), sep = ",")
} else if (group2use=="LANG"){
  lang_eeg_pheno_data = eeg_pheno_data %>% filter(is.element(subid, lang_subs$subid)) %>% distinct()
  write(as.character(lang_eeg_pheno_data$subid), file = file.path(tmp_dir,sprintf("lang_sublist_%s.csv",fstem)), sep = ",")
}
