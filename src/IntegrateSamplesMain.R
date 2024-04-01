#!/usr/local/bin/Rscript
#source(paste0(dirname(dirname(dirname(getwd()))),"/utils/seurat_analysis.R"))

library(argparse)
library(Seurat)
library(dplyr)
library(tidyverse)

set.seed(333)

parser <- ArgumentParser(description='Process scRNA-seq data, while performing comparitive analyses')
parser$add_argument('-input', '--i',  type="character", required=TRUE, nargs=1, help='Contains CellRanger count outputs folder seperated by condition')
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
source(paste0(file.path(dirname(dirname(thisFile()))), "/utils/seurat_analysis.R"))
source(paste0(file.path(dirname(dirname(thisFile()))), "/utils/misc.R"))

se.integrated <- IntegrateSamples(se.filtered.singlets.list, group)
saveRDS(se.integrated, "se_integrated.rds")
