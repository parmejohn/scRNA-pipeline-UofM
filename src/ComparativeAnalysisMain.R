#!/usr/local/bin/Rscript
#source("../utils/deseq2_gsea_analysis.R", chdir=TRUE)

library(argparse)
library(Seurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(ggrepel)
library(DESeq2)
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
))), "/utils/deseq2_gsea_analysis.R"))
source(paste0(file.path(dirname(dirname(
  thisFile()
))), "/utils/misc.R"))


species <- NA
if (args$s == "musmusculus") {
  species <- "Mus musculus"
} else if (args$s == "homosapiens") {
  species <- "Homo sapiens"
} else {
  print("bad")
}

print("DESEQ2 and GSEA")
if (length(unique(se.integrated@meta.data[["group"]])) >= 2) {
  DESeq2ConditionPerCluster(se.integrated, species, args$plots_format)
}
