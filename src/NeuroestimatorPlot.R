#!/usr/local/bin/Rscript

library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(tidyverse)

set.seed(333)
## set command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

#stop the script if no command line argument
if(length(args)==0){
  print("Please include neuroestimator results")
  stop("Requires command line argument.")
}

thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
source(paste0(file.path(dirname(dirname(thisFile()))), 
              "/utils/neuroestimator_plot.R"))
source(paste0(file.path(dirname(dirname(thisFile()))), 
              "/utils/misc.R"))


# args[1]=neuroestimator results
# args[2]=sample table
# args[3]=plots format
NeuroestimatorPlot(args[1], args[2], args[3])
