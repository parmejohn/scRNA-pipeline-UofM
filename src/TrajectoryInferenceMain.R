#!/usr/local/bin/Rscript
source(paste0(dirname(dirname(dirname(getwd()))),"/utils/trajectory_inference_analysis.R"))

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

se.integrated <- readRDS(indir)

se.integrated <- readRDS("/home/projects/sc_pipelines/test_run_nf_1/analysis/data/se_integrated_auto_label.rds")

print("Trajectory Inference")
# https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html#basic_processing_with_seurat_pipeline
if (args$beginning_cluster == '') {
  TrajectoryInferenceSlingshot(se.integrated)
} else {
  TrajectoryInferenceSlingshot(se.integrated, start.clus=args$beginning_cluster)
}
