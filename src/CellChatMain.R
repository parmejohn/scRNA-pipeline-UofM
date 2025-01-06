#!/usr/local/bin/Rscript
#source(paste0(dirname(dirname(dirname(getwd()))),"/utils/trajectory_inference_analysis.R"))

library(argparse)
library(CellChat)
library(Seurat)
library(ComplexHeatmap)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(svglite)

set.seed(333)

parser <- ArgumentParser(description='Process scRNA-seq data, while performing comparitive analyses')
parser$add_argument('-input', '--i',  type="character", required=TRUE, nargs=1, help='Contains CellRanger count outputs folder seperated by condition')
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
source(paste0(file.path(dirname(dirname(thisFile()))), "/utils/cellchat_analysis.R"))
source(paste0(file.path(dirname(dirname(thisFile()))), "/utils/misc.R"))

se.integrated <- readRDS(indir)

print("CellChat")

species <- NA
if (args$s == "musmusculus") {
  species <- "Mus musculus"
} else if (args$s == "homosapiens") {
  species <- "Homo sapiens"
} else {
  print("bad")
}

if (length(unique(se.integrated@meta.data[["group"]])) >= 2) {
  CellChatAnalysis(se.integrated, species, args$plots_format)
}
