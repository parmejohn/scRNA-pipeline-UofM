#!/usr/local/bin/Rscript
#source(paste0(dirname(dirname(dirname(getwd()))),"/utils/seurat_analysis.R"))

library(argparse)
library(Seurat)
library(dplyr)
library(tidyverse)
library(NbClust)

set.seed(333)

parser <-
  ArgumentParser(description = 'Process scRNA-seq data, while performing comparitive analyses')
parser$add_argument(
  '-input',
  '--i',
  type = "character",
  required = TRUE,
  nargs = 1,
  help = 'Contains CellRanger count outputs folder seperated by condition'
)
parser$add_argument('-clusters_optimal',
                    type = "integer",
                    nargs = 1,
                    help = 'Optimal clusters for dimensional reductions and clustering algorithms')
parser$add_argument('-resolution',
                    type = "integer",
                    nargs = 1,
                    help = 'Change resolution for Seurat UMAP, higher values will lead to more clustering')
args <- parser$parse_args()

input <- args$i
se.integrated <- readRDS(input)
opt.clusters <- 0

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
))), "/utils/seurat_analysis.R"))
source(paste0(file.path(dirname(dirname(
  thisFile()
))), "/utils/misc.R"))

print("Dimensional reduction")
if (args$clusters_optimal == 0) {
  opt.clusters <-
    NbClust(
      se.integrated@reductions[["pca"]]@feature.loadings,
      distance = "euclidean",
      min.nc = 10,
      max.nc = 20,
      method = "complete",
      index = "ch"
    )$Best.nc[1] %>% unname()
} else {
  opt.clusters <- args$clusters_optimal
}

if (args$resolution == 1) {
  se.integrated <-
    SeuratDimReduction(se.integrated, 1:opt.clusters, 'group')
} else {
  se.integrated <-
    SeuratDimReduction(se.integrated, 1:opt.clusters, 'group', args$resolution)
}

# have to use visual check to find optimal # of clusters for now
# this can always be increased if want to seperate clusters more ie. finding rare cell populations
elbow <- ElbowPlot(se.integrated) +
  ggtitle("Elbow plot for Determining Optimal Clusters")
PrintSave(elbow, "integrated_elbow_plot.pdf")

#read.table("/home/projects/sc_pipelines/test_run_nf_1/analysis/data/optimal_clusters.txt")
write.table(
  opt.clusters,
  "optimal_clusters.txt",
  quote = FALSE,
  sep = '',
  row.names = F,
  col.names = F,
  eol = ''
)
saveRDS(se.integrated, "se_integrated_dimred.rds")
