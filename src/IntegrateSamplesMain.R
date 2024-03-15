#!/usr/local/bin/Rscript
source(paste0(dirname(dirname(dirname(getwd()))),"/utils/seurat_analysis.R"))

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

se.integrated <- IntegrateSamples(se.filtered.singlets.list, group)
saveRDS(se.integrated, "se_integrated.rds")
