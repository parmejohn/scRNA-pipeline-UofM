#!/usr/local/bin/Rscript
#source(paste0(dirname(dirname(dirname(getwd()))),"/utils/trajectory_inference_analysis.R"))

library(argparse)
library(Seurat)
library(RColorBrewer)
library(slingshot)
library(SingleCellExperiment)
library(tidyverse)
library(tibble)
library(gam)
library(tradeSeq)
library(ComplexHeatmap)
library(magrittr)

set.seed(333)

parser <- ArgumentParser(description='Process scRNA-seq data, while performing comparitive analyses')
parser$add_argument('-inputr', '--i',  type="character", required=TRUE, nargs=1, help='Contains CellRanger count outputs folder seperated by condition')
parser$add_argument('-beginning_cluster', type="character", nargs=1, help='Beginning cluster for trajectory inference with slingshot')
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
source(paste0(file.path(dirname(dirname(thisFile()))), "/utils/trajectory_inference_analysis.R"))
source(paste0(file.path(dirname(dirname(thisFile()))), "/utils/misc.R"))

se.integrated <- readRDS(indir)

print("Trajectory Inference")
# https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html#basic_processing_with_seurat_pipeline
if (args$beginning_cluster == 'none') {
  TrajectoryInferenceSlingshot(se.integrated)
} else {
  TrajectoryInferenceSlingshot(se.integrated, start.clus=args$beginning_cluster)
}
