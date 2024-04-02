#!/usr/local/bin/Rscript
#source("/utils/da_analysis.R", chdir=TRUE)

library(argparse)
library(Seurat)
library(presto)
library(miloR)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(msigdbr)
library(ggrepel)
library(SingleCellExperiment)
library(patchwork)
library(ggbeeswarm)

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
parser$add_argument('-clusters_optimal',
                    type = "integer",
                    nargs = 1,
                    help = 'Optimal clusters for dimensional reductions and clustering algorithms')
args <- parser$parse_args()

indir <- args$i
se.integrated <- readRDS(indir)

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
))), "/utils/da_analysis.R"))
source(paste0(file.path(dirname(dirname(
  thisFile()
))), "/utils/misc.R"))


##### Differential Abundance #####
print("Differential Abundance")

# make sure idents are set correctly
# Idents(object = se.integrated) <- "celltype"
DifferentialAbundanceMilo(
  se.integrated,
  'sample',
  'group',
  k = args$clusters_optimal,
  d = args$clusters_optimal,
  'HARMONY'
)
