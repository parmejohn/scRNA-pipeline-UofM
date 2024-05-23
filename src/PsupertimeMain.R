#!/usr/local/bin/Rscript

library(argparse)
library(psupertime)
library(RColorBrewer)
library(gtools)
library(Matrix)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)

set.seed(333)

parser <- ArgumentParser(description='Process scRNA-seq data, while performing comparitive analyses')
parser$add_argument('-input', '--i',  type="character", required=TRUE, nargs=1, help='Contains CellRanger count outputs folder seperated by condition')
parser$add_argument(
  '-main_time',
  type = "character",
  required = TRUE,
  nargs = 1,
  help = 'Is the main condition time?'
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

source(paste0(file.path(dirname(dirname(thisFile()))), "/utils/psupertime.R"))
source(paste0(file.path(dirname(dirname(thisFile()))), "/utils/misc.R"))

se.integrated <- readRDS(indir)

RunPsuper(se.integrated, args$main_time)