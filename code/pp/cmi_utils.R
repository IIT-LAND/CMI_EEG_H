# cmi_utils.R

#-------------------------------------------------------------------------------
library(ggpackets)

geom_scatterbox <- ggpacket() +
  geom_jitter() +
  geom_boxplot(fill = NA, colour = "#000000", outlier.shape = NA)



#-------------------------------------------------------------------------------
# function to scale variables before modeling
scale_variables <- function(data2use, vars2scale){
  for (curr_var in vars2scale){
    data2use[,curr_var] = scale(data2use[,curr_var])
  }
  return(data2use)
} # scale_variables



#-------------------------------------------------------------------------------
# function to find subjects with a specific kind of diagnosis
dx_finder <- function(pheno_data, dx2find){

  print("selecting columns to use")
  dx_data1 = pheno_data %>% select(Basic_Demos.EID, contains("Diagnosis_ClinicianConsensus"))
  dx_data2 = pheno_data %>% select(Basic_Demos.EID, contains("ConsensusDx"))
  dx_data = merge(dx_data1,dx_data2, by = "Basic_Demos.EID")

  print("pivoting data frame")
  dx_data_long = dx_data %>%
    pivot_longer(c(Diagnosis_ClinicianConsensus.DX_01,
                   Diagnosis_ClinicianConsensus.DX_02,
                   Diagnosis_ClinicianConsensus.DX_03,
                   Diagnosis_ClinicianConsensus.DX_04,
                   Diagnosis_ClinicianConsensus.DX_05,
                   Diagnosis_ClinicianConsensus.DX_06,
                   Diagnosis_ClinicianConsensus.DX_07,
                   Diagnosis_ClinicianConsensus.DX_08,
                   Diagnosis_ClinicianConsensus.DX_09,
                   Diagnosis_ClinicianConsensus.DX_10,
                   ConsensusDx.DX_01,
                   ConsensusDx.DX_02,
                   ConsensusDx.DX_03,
                   ConsensusDx.DX_04,
                   ConsensusDx.DX_05,
                   ConsensusDx.DX_06,
                   ConsensusDx.DX_07,
                   ConsensusDx.DX_08,
                   ConsensusDx.DX_09,
                   ConsensusDx.DX_10)) %>%
    select(Basic_Demos.EID,name,value)

  sublist = unique(dx_data_long$Basic_Demos.EID)
  subs2use = character()
  dx_found = character()
  for (sub in sublist){
    # print(sub)
    sub_data = dx_data_long %>% filter(Basic_Demos.EID==sub)
    diagnoses = unique(sub_data$value)
    dx_test = grep(dx2find, diagnoses, value=TRUE)
    # if (is.element(diagnoses,dx2find)){
    if (!is_empty(dx_test)){
      subs2use = c(subs2use,rep(sub,length(dx_test)))
      dx_found = c(dx_found, dx_test)
    }
  }

  res_df = data.frame(subid = subs2use, dx = dx_found)
  # return(subs2use)
  return(res_df)

} # function dx_finder
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# function to compute effect sizes
cohens_d <- function(x, y, DIM=1, SIGN=TRUE, na.rm=TRUE) {
  #
  # Function will compute cohen's d effect size.
  # Generalized to work on either matrices, data frames, vectors
  #
  # INPUT
  #	x <- matrix or numeric vector
  #	y <- matrix or numeric vector
  #	DIM <- specify the dimension which samples are along
  #
  # Example usage
  #
  # x <- cbind(rnorm(100,0,1), rnorm(100,0,1), rnorm(100,0,1))
  # y <- cbind(rnorm(100,1,1), rnorm(100,2,1), rnorm(100,3,1))
  # d <- cohens_d(x, y, 1)
  #
  # written by mvlombardo - 28.08.2015
  #

  library(matrixStats)

  # if x and y are vectors, coerce them into matrices
  if (class(x)=="numeric" | class(x)=="integer") {
    x <- as.matrix(x)
  } # if

  if (class(y)=="numeric" | class(y)=="integer") {
    y <- as.matrix(y)
  }# if

  if (na.rm==TRUE){
    missingValDecision = TRUE
  } else {
    missingValDecision = FALSE
  }

  # n-1 for x and y
  lx <- dim(x)[DIM]-1
  ly <- dim(y)[DIM]-1

  # if samples are along the rows
  if (DIM==1){
    if (SIGN){
      # mean difference (numerator)
      md <- colMeans(x, na.rm = missingValDecision) - colMeans(y, na.rm = missingValDecision)
    } else{
      # mean difference (numerator)
      md <- abs(colMeans(x, na.rm = missingValDecision) - colMeans(y, na.rm = missingValDecision))
    }# if (SIGN)
    # pooled variance (denominator), but before any sqrt is done
    csd <- (lx * rowVars(t(x),na.rm = missingValDecision)) + (ly * rowVars(t(y), na.rm = missingValDecision))

    # else if samples are along the columns
  } else if (DIM==2){
    if (SIGN){
      # mean difference (numerator)
      md <- rowMeans(x, na.rm = missingValDecision) - rowMeans(y, na.rm = missingValDecision)
    } else{
      # mean difference (numerator)
      md <- abs(rowMeans(x, na.rm = missingValDecision) - rowMeans(y, na.rm = missingValDecision))
    }# if (SIGN)
    # pooled variance (denominator), but before any sqrt is done
    csd <- lx * rowVars(x, na.rm = missingValDecision) + ly * rowVars(y, na.rm = missingValDecision)
  }# end if

  # divide pooled variance by sum of n-1 for x and y and then square root it
  csd <- sqrt(csd/(lx + ly))
  # compute cohen's d
  cd  <- md/csd
}# end cohens_d <- function(x, y, DIM)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# function to get a set of colors for ggplot
# get_ggColorHue.R
#
# Function will spit out color codes for a given n, just like ggplot2.
# Helpful when you want to manually set colors in plots.
#

get_ggColorHue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
} # function get_ggColorHue
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# function for making jitter boxplot
jitterbox <- function(data, x, y, colour, facet=NULL){
  require(tidyverse)

  p = data %>%
    ggplot(aes_string(x = x, y = y, colour=colour))

  if (!is.null(facet)){
    # p = p + facet_grid(cols = vars(facet))
    p = p + facet_wrap(facet)
  }

  p = p + geom_jitter() +
    geom_boxplot(fill = NA, colour = "#000000", outlier.shape = NA)

  return(p)
} # function jitterbox
#-------------------------------------------------------------------------------
