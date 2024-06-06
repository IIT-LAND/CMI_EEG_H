# cmi_tidy_hurst.R

# load packages
library(easypackages)
libraries("here","tidyverse")

tidy_hurst <- function(masterfile, datafile){
  
  # paths
  codepath = here("code")
  datapath = here("data")
  postprocpath = here("data","postproc")
  tidypath = here("data","tidy")
  
  # parse apart datafile to get sub_name and task_name
  metric_path = dirname(datafile)
  task_path = dirname(metric_path)
  metric_name = basename(metric_path)
  task_name = basename(task_path)
  sub_path = dirname(task_path)
  sub_name = basename(sub_path)
  
  # read in datafile
  print(sprintf("Reading in %s", datafile))
  datafile_df = read.csv(file=datafile)
  
  # change second column name to the metric_name, if its not RestingState
  if (task_name!="RestingState"){
    
    colnames(datafile_df)[is.element(colnames(datafile_df),task_name)] = metric_name
    colnames_datafile = colnames(datafile_df)
    
    # add in columns for subid and task
    datafile_df$subid = sub_name
    datafile_df$task = task_name
    
    # rearrange columns so that subid and task are the first two columns
    datafile_df = datafile_df[,c("subid","task",colnames_datafile)]
    
  } else {
    
    datafile_df_long = datafile_df %>% 
      pivot_longer(c("open1","open2","open3","open4","open5",
                     "closed1","closed2","closed3","closed4","closed5"))
    
    # make condition
    open_mask = is.element(datafile_df_long$name,c("open1","open2","open3","open4","open5"))
    closed_mask = is.element(datafile_df_long$name,c("closed1","closed2","closed3","closed4","closed5"))
    datafile_df_long$condition = NA
    datafile_df_long$condition[open_mask] = "open"
    datafile_df_long$condition[closed_mask] = "closed"
    
    # make block
    datafile_df_long$block = NA
    datafile_df_long$block[is.element(datafile_df_long$name,c("open1","closed1"))] = 1
    datafile_df_long$block[is.element(datafile_df_long$name,c("open2","closed2"))] = 2
    datafile_df_long$block[is.element(datafile_df_long$name,c("open3","closed3"))] = 3
    datafile_df_long$block[is.element(datafile_df_long$name,c("open4","closed4"))] = 4
    datafile_df_long$block[is.element(datafile_df_long$name,c("open5","closed5"))] = 5
    
    # change value and name column to metric_name and condition_name
    colnames(datafile_df_long)[colnames(datafile_df_long)=="value"] = metric_name
    colnames(datafile_df_long)[colnames(datafile_df_long)=="name"] = "condition_name"
    
    # grab existing set of column names
    colnames_datafile = colnames(datafile_df_long)
    
    # make subid and task
    datafile_df_long$subid = sub_name
    datafile_df_long$task = task_name
    
    # rearrange columns so that subid and task are the first two columns
    datafile_df_long = datafile_df_long[,c("subid","task",colnames_datafile)]
    
    # replace datafile_df with the new datafile_df_long
    datafile_df = datafile_df_long
    
  } # if (task_name!="RestingState"){
  
  # check if masterfile already exists
  master_file_exists = file.exists(masterfile)
  
  print(sprintf("Adding %s %s data to %s", sub_name, task_name, masterfile))
  if (master_file_exists){
    
    # read in the masterfile
    result_df = read.csv(masterfile, row.name=1)
    
    # do rbind to concatenate with result_df
    result_df = rbind(result_df, datafile_df)
    
  } else {
    
    # result_df is datafile_df
    result_df = datafile_df
    
  } # if (master_file_exists){
  
  # write out the file to disk
  write.csv(result_df, file = masterfile)
  
  return(result_df)
  
} # function tidy_hurst


cmi_tidy_hurst <- function(sublist,task,masterfile){

  sublist = sublist$V1
  for (sub in sublist){
    datafile = here("data","postproc",sub,task,"H",sprintf("%s_%s_H.csv",sub,task))
    result_df = tidy_hurst(masterfile=masterfile, datafile=datafile)
  } # for (sub in sublist){
  
  return(result_df) 
  
} # cmi_tidy_hurst


