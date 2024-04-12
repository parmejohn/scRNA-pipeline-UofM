#' ---
#' title: "scRNA-seq analysis pipeline report"
#' params:
#'   data:
#'   opt_c:
#' output:
#'   html_document:
#'    toc: true
#'    number_sections: true
#'    toc_float: yes
#'    code_folding: show
#' ---
#'
#'
#' <style type="text/css">
#'   .main-container {
#'     max-width: 100% !important;
#'     margin: auto;
#'   }
#' </style>
  

#+ echo=FALSE, include=FALSE
#library(argparse)
library(dplyr)
library(tidyverse)
library(Seurat)
library(knitr)
library(data.table)
#library(kableExtra)
library(cowplot)
library(magick)

set.seed(333)

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

# use Rmd params option for rendering https://rmarkdown.rstudio.com/lesson-6.html
#res.loc <- "/home/projects/sc_pipelines/scrna_deanne_harmony/downsample_test/analysis/"
res.loc <- params$data

#' # QC
#' ## Ambient RNA Removal
#' List and location of ambient RNA QC'd files
#+ warning=FALSE, echo=FALSE
dir1 <- list.dirs(paste0(res.loc, "/data/qc"), recursive = FALSE)
dir2 <- list.dirs(dir1, recursive = FALSE)
print(dir2)

#' ## Filtering Low Quality Cells
#' Doublets are also removed, but no accompanying graph follows. 
#+ warning=FALSE, echo=FALSE, fig.show='hold', fig.height = 10, fig.width = 10, out.width="50%", results='asis'
plots <- list.files(paste0(res.loc, "plots/qc/"), full.names=TRUE)

for(i in plots){
  #cat("![example](",i,"){width=100%, height=500}")
  #cat("\includepdf[pages={1}](",i,"){width=100%}")
  #cat("<iframe src=",i,"width='100%' height='500' frameborder='0' />")
  plot(image_read_pdf(i))
}

#' # Seurat processing
#' ## Integrating and Performing dimensional reductions
#' Note that this integration will be done on the following condtions
{{list.dirs(paste0(res.loc, "/data/qc"), recursive = FALSE, full.names = FALSE)}}
#' , which will be named 'group' in the seurat metadata
#' Optimal clustering was inputed as 
{{params$optc}}

#' Please view the grouped plot below to ensure that integration was sufficient for your given dataset
#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10

plot(image_read_pdf(paste0(res.loc,"plots/integrated_umap_unlabelled.pdf")))
plot(image_read_pdf(paste0(res.loc,"plots/integrated_umap_grouped.pdf")))
plot(image_read_pdf(paste0(res.loc,"plots/integrated_umap_split.pdf")))

se.integrated <- readRDS(paste0(res.loc, "data/se_integrated_dimred.rds"))
md <- se.integrated@meta.data %>% as.data.table()

kable(md[, .N, by = c("group", "seurat_clusters")] %>% dcast(., group ~ seurat_clusters, value.var = "N"), caption = "Number of cells per seurat_cluster")

#' ## Identifying Cell Markers
#' Please use the information below to identify/confirm your clustering labels.
#' The heatmap below is also limited to 300 cells per cluster because of ggplot limitations.
#' Note: Cell marker in this context means genes that are differentially expressed in the given cluster compared to all others.
#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10, results='asis'
plot(image_read_pdf(paste0(res.loc,"plots/top3_markers_5k_expr_heatmap.pdf")))

if (file.exists(paste0(res.loc,"plots/conserved_marker_unlabelled.pdf"))){
  plot(image_read_pdf(paste0(res.loc,"plots/conserved_marker_unlabelled.pdf")))
}

kable(head(read.table(paste0(res.loc,"data/se_markers_presto_integrated.txt"))), caption="Top cell markers for each cluster (showing just 10)")

#' ### Automatic Clustering Labelling
#' Remember to double check how accurate the automatic labelling was!
#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10,results='asis'
if (file.exists(paste0(res.loc,"data/se_integrated_auto_label.rds"))){

  plot(image_read_pdf(paste0(res.loc,"plots/reference_marker_mapping_heatmap.pdf")))
  plot(image_read_pdf(paste0(res.loc,"plots/integrated_umap_labelled.pdf")))
  
  se.integrated <- readRDS(paste0(res.loc, "data/se_integrated_auto_label.rds"))
  
  se.integrated$labelled <- Idents(se.integrated)
  
  md <- se.integrated@meta.data %>% as.data.table()
  kable(md[, .N, by = c("group", "labelled")] %>% dcast(., group ~ labelled, value.var = "N"), caption = "Number of cells per labelled cluster")
} else {
  print("No reference seurat was specified, rest of the analysis will use the default cluster labelling from seurat")
}

#' # Post-Clustering Analysis
#' ## Comparative Analysis
#' The following graphs use an aggregated expression from each cluster for comparisons.
#' 
#' ### Differential Gene Expression Analysis with DESeq2
#' Please note that only 10 comparisons were chosen and printed here, if you would like to see all the DESeq2 plots generated, go to
{{paste0(res.loc, "plots/deseq2")}}

#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10, out.width = "50%", fig.show = 'hold', results='asis'
plots <- list.files(paste0(res.loc, "plots/deseq2/"), full.names=TRUE)

for(i in 1:length(plots)){
  #cat("![example](",i,"){width=100%, height=500}")
  #cat("\includepdf[pages={1}](",i,"){width=100%}")
  if (i != 10){
    plot(image_read_pdf(plots[i]))
  } else {
    break
  }
}

#' ### GSEA Analysis
#' This uses the DESeq2 results for the ranked gene list, which is used to calculate the enrichment score
#' All the GSEA plots are available at 
{{paste0(res.loc, "plots/gsea/comparative")}}

#+ warning=FALSE, echo=FALSE, fig.height = 8, fig.width = 10, out.width = "50%", fig.show = 'hold', results='asis'
plots <- list.files(paste0(res.loc, "plots/gsea/comparative/"), full.names=TRUE)

for(i in 1:length(plots)){
  #cat("![example](",i,"){width=100%, height=500}")
  #cat("\includepdf[pages={1}](",i,"){width=100%}")
  if (i != 10){
    plot(image_read_pdf(plots[i]))
  } else {
    break
  }
}

#' ## Single-cell GSEA with escape
#' Printed is the top 5 average normalized enrichment scores from the escape analysis.
#' Individual GO pathways in a geyser plot can be found at 
{{paste0(res.loc, "plots/gsea/escape")}}

#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10, results='asis'
if (file.exists(paste0(res.loc,"data/se_integrated_escape_norm.rds"))){
  plot(image_read_pdf(paste0(res.loc,"plots/gsea/escape/escape_heatmap_top5.pdf")))
} else {
  print("escape analysis was set to false, if you wanted this analysis please set '-run_escape true'")
}

#' ## Trajectory Inference with slingshot
#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10, results='asis'
if (file.exists(paste0(res.loc,"plots/ti"))){
  if(file.exists(paste0(res.loc, "plots/ti/ti_start_smooth.pdf"))){
    plots <- list.files(paste0(res.loc, "plots/ti/"), full.names=TRUE)
    for(i in plots){
      #cat("![example](",i,"){width=100%, height=500}")
      #cat("\includepdf[pages={1}](",i,"){width=100%}")
#      plot(image_read_pdf(i))
    }
  } else {
#    plot(image_read_pdf(paste0(res.loc,"plots/ti/ti_no_start_not_smooth.pdf")))
  }
} else {
  print("slingshot analysis was set to false, if you wanted this analysis please set '-run_sling true'")
}

#' ## Differential Abundance analysis with miloR
{{if(length(list.files(paste0(res.loc, "plots/da/"), full.names=TRUE)) == 3){"No Differential Abundant neighborhoods were found, so no DEG heatmaps will be made."}}}
#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10, results='asis'
plots <- list.files(paste0(res.loc, "plots/da/"), full.names=TRUE)

for(i in plots){
  #cat("![example](",i,"){width=100%, height=500}")
  #cat("\includepdf[pages={1}](",i,"){width=100%}")
  plot(image_read_pdf(i))
}

#knitr::spin(thisFile())
#knitr::spin("/home/projects/sc_pipelines/scRNA-pipeline-UofM_harmony/src/SummaryReport.R")
# rmarkdown::render("/home/projects/sc_pipelines/scRNA-pipeline-UofM_harmony/src/SummaryReport.R")
