#!/usr/local/bin/Rscript
#source(paste0(dirname(dirname(dirname(getwd()))),"/utils/cellmarker_analysis.R"))

library(argparse)
library(Seurat)
library(SeuratWrappers)
library(presto)
library(dplyr)
library(tidyverse)
library(pheatmap)

set.seed(333)

parser <- ArgumentParser(description='Process scRNA-seq data, while performing comparitive analyses')
parser$add_argument('-input', '--i',  type="character", required=TRUE, nargs=1, help='Contains CellRanger count outputs folder seperated by condition')
parser$add_argument('-clusters_optimal', type="integer", nargs=1, help='Optimal clusters for dimensional reductions and clustering algorithms')
parser$add_argument('-reference_seurat', type="character", nargs=1, help='Reference Seurat(s) object file path with pre-labelled clusters') # need to change to many options later
args <- parser$parse_args()

input <- args$i
se.integrated <- readRDS(input)

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
source(paste0(file.path(dirname(dirname(thisFile()))), "/utils/cellmarker_analysis.R"))
source(paste0(file.path(dirname(dirname(thisFile()))), "/utils/misc.R"))

#### Identifying Cell Markers ####
print("IdentifyCellMarkers")
IdentifyCellMarkers(se.integrated)

##### external dataset #####
print("External mapping")
if (args$reference_seurat != "none") {
  se.reference <- readRDS(args$reference_seurat)
  se.integrated <- ReferenceMarkerMapping(se.reference, se.integrated, args$clusters_optimal) # active.ident should correlate to the cluster identification
  saveRDS(se.integrated, "se_integrated_auto_label.rds")
}
