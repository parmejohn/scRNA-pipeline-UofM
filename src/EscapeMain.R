#!/usr/local/bin/Rscript
#source("../utils/sc_gsea.R", chdir=TRUE)

library(argparse)
library(Seurat)
library(presto)
library(escape)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(msigdbr)
library(ggrepel)

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
))), "/utils/sc_gsea.R"))
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

print("single-cell GSEA")
EscapeGSEA(se.integrated, species)
