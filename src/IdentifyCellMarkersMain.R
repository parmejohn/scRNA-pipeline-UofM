#!/usr/local/bin/Rscript
source(paste0(dirname(dirname(dirname(getwd()))),"/utils/cellmarker_analysis.R"))

library(argparse)
library(Seurat)
library(dplyr)
library(tidyverse)

set.seed(333)

parser <- ArgumentParser(description='Process scRNA-seq data, while performing comparitive analyses')
parser$add_argument('-input', '--i',  type="character", required=TRUE, nargs=1, help='Contains CellRanger count outputs folder seperated by condition')
parser$add_argument('-clusters_optimal', type="integer", nargs=1, help='Optimal clusters for dimensional reductions and clustering algorithms')
parser$add_argument('-resolution', type="integer", nargs=1, help='Change resolution for Seurat UMAP, higher values will lead to more clustering')
args <- parser$parse_args()

input <- args$i