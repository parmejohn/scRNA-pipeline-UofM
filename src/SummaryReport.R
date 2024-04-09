#' ---
#' title: "scRNA-seq analysis pipeline report"
#' output:
#'   html_document:
#'    toc: true
#'    number_sections: true
#'    toc_float: yes
#'    code_folding: show
#' ---
#'
#'
#+ echo=FALSE, include=FALSE
library(argparse)
library(dplyr)
library(tidyverse)
library(Seurat)
library(knitr)
library(data.table)
#library(kableExtra)

#library(magick)

set.seed(333)

parser <-
  ArgumentParser(description = 'Process scRNA-seq data, while performing comparitive analyses')
parser$add_argument(
  '-analysis_folder',
  '--a',
  type = "character",
  required = TRUE,
  nargs = 1,
  help = 'Contains pipeline analysis folder path'
)

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

res.loc <- "/home/projects/sc_pipelines/scrna_deanne_harmony/downsample_test/analysis/"

#' # QC
#' ## Ambient RNA Removal
#' List and location of ambient RNA QC'd files
#+ warning=FALSE, echo=FALSE
dir1 <- list.dirs(paste0(res.loc, "/data/qc"), recursive = FALSE)
dir2 <- list.dirs(dir1, recursive = FALSE)
print(dir2)

#' ## Filtering Low Quality Cells
#' Doublets are also removed, but no accompanying graph follows. 
#+ warning=FALSE, echo=FALSE, results='asis'
plots <- list.files(paste0(res.loc, "plots/qc/"), full.names=TRUE)

for(i in plots){
  #cat("![example](",i,"){width=100%, height=500}")
  #cat("\includepdf[pages={1}](",i,"){width=100%}")
  cat("<iframe src=",i,"width='100%' height='500' frameborder='0' />")
}

#' # Seurat processing
#' ## Integrating and Performing dimensional reductions
#' Note that this integration will be done on the following condtions, which will be named 'group' in the seurat metadata
{{list.dirs(paste0(res.loc, "/data/qc"), recursive = FALSE, full.names = FALSE)}}
#' Optimal clustering was inputed as 
#{{args$optimal_clusters}}

#' Please view the grouped plot below to ensure that integration was sufficient for your given dataset
#+ warning=FALSE, echo=FALSE, results='asis'
cat("![example](",res.loc,"plots/integrated_umap_unlabelled.pdf){width=100%, height=500}", sep="")
cat("![example](",res.loc,"plots/integrated_umap_grouped.pdf){width=100%, height=500}", sep="")
cat("![example](",res.loc,"plots/integrated_umap_split.pdf){width=100%, height=500}", sep="")

se.integrated <- readRDS(paste0(res.loc, "data/se_integrated_dimred.rds"))
md <- se.integrated@meta.data %>% as.data.table()

kable(md[, .N, by = c("group", "seurat_clusters")] %>% dcast(., group ~ seurat_clusters, value.var = "N"), caption = "Number of cells per seurat_cluster")

#' ## Identifying Cell Markers
#' Please use the information below to identify/confirm your clustering labels.
#' The heatmap below is also limited to
#+ warning=FALSE, echo=FALSE, results='asis'
cat("![example](",res.loc,"plots/top3_markers_expr_heatmap.pdf){width=100%, height=500}", sep="")

if (file.exists(paste0(res.loc,"plots/conserved_marker_unlabelled.pdf"))){
  cat("![example](",res.loc,"plots/conserved_marker_unlabelled.pdf){width=100%, height=500}", sep="")
}

kable(head(read.table(paste0(res.loc,"data/se_markers_presto_integrated.txt"))), caption="Top cell markers for each cluster (showing just 10)")

#' ### Automatic Clustering Labelling
#' Remember to double check how accurate the automatic labelling was!
#+ warning=FALSE, echo=FALSE, results='asis'
if (file.exists(paste0(res.loc,"data/se_integrated_auto_label.rds"))){
  cat("![example](",res.loc,"plots/integrated_umap_labelled.pdf){width=100%, height=500}", sep="")
  se.integrated <- readRDS(paste0(res.loc, "data/se_integrated_auto_label.rds"))
  
  se.integrated$labelled <- Idents(se.integrated)
  
  md <- se.integrated@meta.data %>% as.data.table()
  kable(md[, .N, by = c("group", "labelled")] %>% dcast(., group ~ labelled, value.var = "N"), caption = "Number of cells per labelled cluster")
} else {
  print("No reference seurat was specified, rest of the analysis will use the default cluster labelling from seurat")
}

#knitr::spin(thisFile())
#knitr::spin("/home/projects/sc_pipelines/scRNA-pipeline-UofM_harmony/src/SummaryReport.R")
#rmarkdown::render("/home/projects/sc_pipelines/scRNA-pipeline-UofM_harmony/src/SummaryReport.R")
