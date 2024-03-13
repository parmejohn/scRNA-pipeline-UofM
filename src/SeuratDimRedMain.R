#!/usr/local/bin/Rscript
source(paste0(dirname(dirname(dirname(getwd()))),"/utils/seurat_analysis.R"))

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

print("Dimensional reduction")
if (length(args$clusters_optimal) == 0) {
  opt.clusters <- NbClust(se.integrated@reductions[["integrated.cca"]]@feature.loadings, distance = "euclidean", min.nc=10, max.nc = 20,
                          method = "complete", index = "ch")$Best.nc[1] %>% unname()
} else {
  opt.clusters <- args$clusters_optimal
}

integrated_elbow <- ElbowPlot(input) # have to use visual check to find optimal # of clusters for now
# this can always be increased if want to seperate clusters more ie. finding rare cell populations
PrintSave(integrated_elbow, "integrated_elbow_plot.pdf")

if (length(args$resolution) == 1) {
  se.integrated <- SeuratDimReduction(se.integrated, 1:opt.clusters, 'group')
} else {
  se.integrated <- SeuratDimReduction(se.integrated, 1:opt.clusters, 'group', args$resolution)
}

saveRDS(se.integrated, "se_integrated_dimred.rds")