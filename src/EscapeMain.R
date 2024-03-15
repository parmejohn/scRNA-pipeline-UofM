#!/usr/local/bin/Rscript
source(paste0(dirname(dirname(dirname(getwd()))),"/utils/sc_gsea.R"))

library(argparse)
library(Seurat)
library(presto)
library(escape)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(msigdbr)
library(ggrepel)

set.seed(333)

parser <- ArgumentParser(description='Process scRNA-seq data, while performing comparitive analyses')
parser$add_argument('-input', '--i',  type="character", required=TRUE, nargs=1, help='Contains CellRanger count outputs folder seperated by condition')
parser$add_argument('-species', '--s',type="character", required=TRUE, nargs=1, help='Species name (Mus musculus, Homo sapiens); CASE-SENSITIVE')
args <- parser$parse_args()

indir <- args$i
se.integrated <- readRDS(indir)

species <- NA
if (args$s == "musmusculus"){
  species <- "Mus musculus"
} else {
  print("bad")
}

print("single-cell GSEA")
EscapeGSEA(se.integrated, species)