#!/usr/local/bin/Rscript
source(paste0(dirname(dirname(dirname(getwd()))),"/utils/da_analysis.R"))

library(argparse)
library(Seurat)
library(presto)
library(miloR)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(msigdbr)
library(ggrepel)
library(SingleCellExperiment)

set.seed(333)

parser <- ArgumentParser(description='Process scRNA-seq data, while performing comparitive analyses')
parser$add_argument('-input', '--i',  type="character", required=TRUE, nargs=1, help='Contains CellRanger count outputs folder seperated by condition')
parser$add_argument('-clusters_optimal', type="integer", nargs=1, help='Optimal clusters for dimensional reductions and clustering algorithms')
args <- parser$parse_args()

indir <- args$i
se.integrated <- readRDS(indir)

##### Differential Abundance #####
print("Differential Abundance")

# make sure idents are set correctly
# Idents(object = se.integrated) <- "celltype"
DifferentialAbundanceMilo(se.integrated, 'sample', 'group', k = args$clusters_optimal, d = 50, 'INTEGRATED.CCA')