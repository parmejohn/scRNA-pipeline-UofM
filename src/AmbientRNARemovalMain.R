#!/usr/local/bin/Rscript
#print(dirname(dirname(dirname(getwd()))))

source(paste0(dirname(dirname(dirname(getwd()))),"/utils/qc.R"))

library(argparse)
library(Seurat)
library(SoupX)
library(DropletUtils)
#library(glmGamPoi)

set.seed(333)

parser <- ArgumentParser(description='Process scRNA-seq data, while performing comparitive analyses')
parser$add_argument('-indir', '--i',  type="character", required=TRUE, nargs=1, help='Contains CellRanger count outputs folder seperated by condition')
args <- parser$parse_args()

indir <- args$i

print("SoupX filtering")
print(indir)
# loading in the files initially; need to test runtime
filenames <- list.files(path = indir, pattern = "matrix.h5", full.names = TRUE, recursive = T, include.dirs = T)
print(filenames)
filenames <- filenames[grepl("outs", filenames)]

list_of_pairs <- list()
for (i in 1:length(filenames)){
  if (i %% 2 == 1){
    temp_list <- list(c(filenames[i], filenames[i+1]))
    list_of_pairs <- append(list_of_pairs, temp_list)
  }
}

lapply(list_of_pairs, AmbientRNARemoval)
