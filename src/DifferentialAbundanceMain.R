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
library(gtools)
library(purrr)
library(rlang)
library(fgsea)
library(msigdbr)
library(gridExtra)
library(ComplexHeatmap)
library(svglite)

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
parser$add_argument('-reduced_dim',
                    type = "character",
                    nargs = 1,
                    help = 'Reduced dimensions used')
parser$add_argument(
  '-species',
  '--s',
  type = "character",
  required = TRUE,
  nargs = 1,
  help = 'Species name (Mus musculus, Homo sapiens); CASE-SENSITIVE'
)
parser$add_argument(
  '-plots_format',
  type = "character",
  required = TRUE,
  nargs = 1
)
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
))), "/utils/miloR_fixes.R"))
source(paste0(file.path(dirname(dirname(
  thisFile()
))), "/utils/da_analysis.R"))
source(paste0(file.path(dirname(dirname(
  thisFile()
))), "/utils/misc.R"))


##### Differential Abundance #####
print("Differential Abundance")

species <- NA
if (args$s == "musmusculus") {
  species <- "Mus musculus"
} else if (args$s == "homosapiens") {
  species <- "Homo sapiens"
} else {
  print("bad")
}

# make sure idents are set correctly
# Idents(object = se.integrated) <- "celltype"
DifferentialAbundanceMilo(
  se.integrated = se.integrated,
  sample = 'sample',
  k = args$clusters_optimal,
  d = 50,
  reduced.dims = toupper(args$reduced_dim),
  species = species,
  plots.format = args$plots_format
)
