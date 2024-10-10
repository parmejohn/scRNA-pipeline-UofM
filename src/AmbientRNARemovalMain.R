#!/usr/local/bin/Rscript
#print(dirname(dirname(dirname(getwd()))))


library(argparse)
library(Seurat)
library(SoupX)
library(DropletUtils)
#library(glmGamPoi)

set.seed(333)

parser <-
  ArgumentParser(description = 'Process scRNA-seq data, while performing comparitive analyses')
parser$add_argument(
  '-indir',
  '--i',
  type = "character",
  required = TRUE,
  nargs = 1,
  help = 'Contains CellRanger count outputs folder seperated by condition'
)
parser$add_argument('-test_data',
                    type = "integer",
                    nargs = 1,
                    help = 'specify whether this is just a test run, if so your data will be downsampled 3000 cells per sample') # need to change to many options later
args <- parser$parse_args()

indir <- args$i

##### load in helper functions #####
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

source(paste0(file.path(dirname(dirname(
  thisFile()
))), "/utils/qc.R"))
source(paste0(file.path(dirname(dirname(
  thisFile()
))), "/utils/misc.R"))

print("SoupX filtering")
print(indir)
# loading in the files initially; need to test runtime
filenames <-
  list.files(
    path = indir,
    pattern = "feature_bc_matrix.h5",
    full.names = TRUE,
    recursive = T,
    include.dirs = T
  )
#print(filenames)
filenames <- filenames[grepl("outs", filenames)]

list_of_pairs <- list()
for (i in 1:length(filenames)) {
  if (i %% 2 == 1) {
    temp_list <- list(c(filenames[i], filenames[i + 1]))
    list_of_pairs <- append(list_of_pairs, temp_list)
  }
}

lapply(list_of_pairs, AmbientRNARemoval, test=args$test_data)
