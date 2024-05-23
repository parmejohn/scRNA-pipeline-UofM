#!/usr/local/bin/Rscript

library(argparse)
library(Tempora)
library(tidyverse)
library(dplyr)
library(gtools)
library(reshape2)
library(Seurat)
library(SeuratObject)
library(RCurl)
library(extrafont)
library(igraph)
library(RColorBrewer)
library(Matrix)
library(BiocParallel)

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
  '-main_time',
  type = "character",
  required = TRUE,
  nargs = 1,
  help = 'Is the main condition time?'
)
args <- parser$parse_args()

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
))), "/utils/tempora.R"))
source(paste0(file.path(dirname(dirname(
  thisFile()
))), "/utils/tempora_fixes.R"))
source(paste0(file.path(dirname(dirname(
  thisFile()
))), "/utils/misc.R"))

input <- args$i
se.integrated <- readRDS(input)

se.integrated.tempora.seurat.v3 <- RunTempora(se.integrated, args$main_time)
saveRDS(se.integrated.tempora.seurat.v3, "se_integrated_tempora_seurat_v3.rds")