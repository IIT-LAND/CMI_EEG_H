# cmi_tidy_dataquality.R

# load packages
library(easypackages)
libraries("here","tidyverse")

tidy_dataquality <- function(masterfile,
                             dataQualFile,
                             scoreEpochsFile,
                             meanTSNRfile,
                             badChanSampFile){

  # paths
  codepath = here("code")
  datapath = here("data")
  preprocpath = here("data","preproc")
  tidypath = here("data","tidy")

  # parse apart datafile to get sub_name and task_name
  task_path = dirname(dataQualFile)
  task_name = basename(task_path)
  sub_path = dirname(task_path)
  sub_name = basename(sub_path)

  # read in files
  print(sprintf("Reading in data quality measures for %s %s", sub_name, task_name))
  dataQual_df = read.csv(file=dataQualFile)
  scoreEpochs_df = read.csv(file=scoreEpochsFile)
  meanTSNR_df = read.csv(file=meanTSNRfile)
  badChanSamp_df = read.csv(file=badChanSampFile)

  dataQual_df = dataQual_df %>% select(contains("preproc"))
  colnames(dataQual_df)[is.element(colnames(dataQual_df),"preproc_OHA")] = "OHA"
  colnames(dataQual_df)[is.element(colnames(dataQual_df),"preproc_THV")] = "THV"
  colnames(dataQual_df)[is.element(colnames(dataQual_df),"preproc_CHV")] = "CHV"
  colnames(dataQual_df)[is.element(colnames(dataQual_df),"preproc_MAV")] = "MAV"
  colnames(scoreEpochs_df)[is.element(colnames(scoreEpochs_df),"Score")] = "avgEpochScore"

  # cbind together data quality, scorepochs, and mean_tSNR
  datafile_df = cbind(badChanSamp_df,dataQual_df, scoreEpochs_df, meanTSNR_df)
  # datafile_df = cbind(dataQual_df, scoreEpochs_df, meanTSNR_df)
  colnames_datafile = colnames(datafile_df)

  # # make subid and task column
  # datafile_df$subid = sub_name
  # datafile_df$task = task_name
  # datafile_df = datafile_df[,c('subid',"task",colnames_datafile)]

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

} # function tidy_dataquality


cmi_tidy_dataquality <- function(sublist,task,masterfile){

  sublist = sublist$V1
  for (sub in sublist){

    dataQualFile = here("data","preproc",sub,task,sprintf("%s_%s_preproc_dataQuality.csv",sub,task))
    scoreEpochsFile = here("data","preproc",sub,task,sprintf("%s_%s_preproc_scorepochs_summaryEpochScore.csv",sub,task))
    meanTSNRfile = here("data","preproc",sub,task,sprintf("%s_%s_preproc_mean_tsnr.csv",sub,task))
    badChanSampFile = here("data","preproc",sub,task,sprintf("%s_%s_preproc_badSamplesChannels.csv",sub,task))

    possibleError = tryCatch({result_df = tidy_dataquality(masterfile=masterfile,
                                                           dataQualFile=dataQualFile,
                                                           scoreEpochsFile=scoreEpochsFile,
                                                           meanTSNRfile=meanTSNRfile,
                                                           badChanSampFile=badChanSampFile)},
                             error = function(e) e)
    class(possibleError)
    if (inherits(possibleError,"error")){
      next
    }

  } # for (sub in sublist){

  return(result_df)

} # cmi_tidy_dataquality
