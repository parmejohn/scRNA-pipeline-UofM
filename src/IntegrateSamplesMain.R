#!/usr/local/bin/Rscript

library(argparse)
library(Seurat)
library(dplyr)
library(tidyverse)
library(harmony)
library(batchelor)
library(SeuratWrappers)
library(Signac)

set.seed(333)

parser <-
  ArgumentParser(description = 'Process scRNA-seq data, while performing comparitive analyses')
parser$add_argument(
  '-input',
  '--i',
  type = "character",
  required = TRUE,
  nargs = 1,
  help = 'Contains CellRanger count outputs folder seperated by condition'
)
parser$add_argument('-reduced_dim',
                    type = "character",
                    nargs = 1,
                    help = 'Reduced dimensions used')
parser$add_argument(
  '-coconditions',
  type = "character",
  required = TRUE,
  nargs = '*',
  help = 'Co-conditions listed'
)
args <- parser$parse_args()

input <- args$i
se.filtered.singlets.list <- readRDS(input)

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
source(paste0(file.path(dirname(dirname(
  thisFile()
))), "/utils/seurat_analysis.R"))
source(paste0(file.path(dirname(dirname(
  thisFile()
))), "/utils/misc.R"))

se.integrated <- IntegrateSamples(se.filtered.singlets.list, 
                                  args$reduced_dim)
if (args$coconditions[1] != 'none'){
  Misc(se.integrated, slot = "co.conditions") <- args$coconditions
}
saveRDS(se.integrated, "se_integrated.rds")
