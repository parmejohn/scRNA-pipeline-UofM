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
library(dplyr)
library(tidyverse)
library(Seurat)
library(knitr)
library(data.table)
#library(kableExtra)
library(cowplot)
library(magick)

set.seed(333)

ReadImageAndTrim <- function(path) {
  image_read_pdf(path) %>% 
    image_trim() %>%
    plot()
}

# use Rmd params option for rendering https://rmarkdown.rstudio.com/lesson-6.html
#res.loc <- "/home/projects/sc_pipelines/organoid_work/organoid_dn_sampled_v3/pipeline/analysis/"
res.loc <- params$data

#' # QC
#' ## Ambient RNA Removal
#' Reminder that this step adjusts the counts accordingly due to removing ambient RNA, it does not remove cells themselves. List and location of ambient RNA QC'd files
#+ warning=FALSE, echo=FALSE
dir1 <- list.dirs(paste0(res.loc, "/data/qc"), recursive = FALSE)
dir2 <- list.dirs(dir1, recursive = FALSE)
print(dir2)

#' ## Filtering Low Quality Cells
#' Doublets are also removed, but no accompanying graph follows. 
#+ warning=FALSE, echo=FALSE, fig.show='hold', fig.height = 10, fig.width = 10, out.width="50%", results='asis'
plots <- list.files(paste0(res.loc, "plots/qc/"), full.names=TRUE)

for(i in plots){
  ReadImageAndTrim(i)
}

doub.table <-  data.frame(row.names = c("singlet", "doublet"))

se.filtered.doublets.list <- readRDS(paste0(res.loc, "/data/se_filtered_doublets_list.rds"))
for (i in 1:length(se.filtered.doublets.list)){
  singlets <- subset(se.filtered.doublets.list[[i]], subset = scDblFinder.class  == "singlet") %>% colnames() %>% length()
  doublets <- subset(se.filtered.doublets.list[[i]], subset = scDblFinder.class  == "doublet") %>% colnames() %>% length()
  
  doub.table[, levels(se.filtered.doublets.list[[i]]@meta.data[["ident"]])] <- c(singlets, doublets)
}

kable(doub.table, caption = "Number of Doublets Detected with scDblFinder")

#' # Seurat processing
#' ## Integrating and Performing dimensional reductions
#' Note that this integration will be done on the following condtions
{{list.dirs(paste0(res.loc, "/data/qc"), recursive = FALSE, full.names = FALSE)}}
#' , which will be named 'group' in the seurat metadata
#' Optimal clustering was inputed as 
{{params$opt_c}}

#' Please view the grouped plot below to ensure that integration was sufficient for your given dataset
#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10
ReadImageAndTrim(paste0(res.loc,"plots/integrated_umap_unlabelled.pdf"))

se.integrated <- readRDS(paste0(res.loc, "data/se_integrated_dimred.rds"))
md <- se.integrated@meta.data %>% as.data.table()

ReadImageAndTrim(paste0(res.loc,"plots/percent_cells_group_unlabelled.pdf"))

kable(md[, .N, by = c("group", "seurat_clusters")] %>% dcast(., group ~ seurat_clusters, value.var = "N"), caption = "Number of cells per seurat_cluster")

ReadImageAndTrim(paste0(res.loc,"plots/integrated_umap_grouped.pdf"))
ReadImageAndTrim(paste0(res.loc,"plots/integrated_umap_split.pdf"))

#' ## Identifying Cell Markers
#' Please use the information below to identify/confirm your clustering labels.
#' The heatmap below is also limited to 300 cells per cluster because of ggplot limitations.
#' Note: Cell marker in this context means genes that are differentially expressed in the given cluster compared to all others.
#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10, results='asis'
ReadImageAndTrim(paste0(res.loc,"plots/top3_markers_expr_heatmap.pdf"))

if (file.exists(paste0(res.loc,"plots/conserved_marker_unlabelled.pdf"))){
  ReadImageAndTrim(paste0(res.loc,"plots/conserved_marker_unlabelled.pdf"))
}

kable(head(read.table(paste0(res.loc,"data/se_markers_presto_integrated.txt"))), caption="Top cell markers for each cluster (showing just 10)")

#' ### Automatic Clustering Labelling
#' Remember to double check how accurate the automatic labelling was!
#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10,results='asis'
if (file.exists(paste0(res.loc,"data/se_integrated_auto_label.rds"))){

  ReadImageAndTrim(paste0(res.loc,"plots/reference_marker_mapping_heatmap.pdf"))
  ReadImageAndTrim(paste0(res.loc,"plots/integrated_umap_labelled.pdf"))

  ReadImageAndTrim(paste0(res.loc,"plots/percent_cells_group_labelled.pdf"))
  
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
#' Cutoffs: p_val_adj < 0.05 & |avg_log2FC| >= 2
#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10, out.width = "50%", fig.show = 'hold', results='asis'
plots <- list.files(paste0(res.loc, "plots/deseq2/"), full.names=TRUE)

for(i in 1:length(plots)){
  if (i <= 10){
    ReadImageAndTrim(plots[i])
  } else {
    break
  }
}

#' ### GSEA Analysis
#' This uses the DESeq2 results for the ranked gene list, which is used to calculate the enrichment score  
#' All the GSEA plots are available at 
{{paste0(res.loc, "plots/gsea/comparative")}}  
#' Cutoffs: p_val_adj < 0.05
#+ warning=FALSE, echo=FALSE, fig.height = 8, fig.width = 10, out.width = "50%", fig.show = 'hold', results='asis'
plots <- list.files(paste0(res.loc, "plots/gsea/comparative/"), full.names=TRUE)

for(i in 1:length(plots)){
  if (i <= 10){
    ReadImageAndTrim(plots[i])
  } else {
    break
  }
}

#' ## Single-cell GSEA with escape
#' Printed is the top 5 average normalized enrichment scores from the escape analysis.  
#' Individual GO pathways in a geyser plot can be found at 
{{paste0(res.loc, "plots/gsea/escape")}}  
#' Cutoffs: p_val_adj <= 0.05
#+ warning=FALSE, echo=FALSE, fig.height = 8, fig.width = 9, out.width = "100%", results='asis'
plots <- list.files(paste0(res.loc, "plots/gsea/escape/"), full.names=TRUE)
plots <- plots[grep(".pdf", plots)]

geyser_plots <- list.files(paste0(res.loc, "plots/gsea/escape/"), full.names=TRUE, recursive = TRUE)
geyser_plots <- geyser_plots[grep("geyser", geyser_plots)]

if (file.exists(paste0(res.loc,"data/se_integrated_escape_norm.rds"))){
  for(i in plots){
    ReadImageAndTrim(i)
  }
  
  for (i in 1:5){
    ReadImageAndTrim(geyser_plots[i])
  }
  
} else {
  print("escape analysis was set to false, if you wanted this analysis please set '-run_escape true'")
}

#' ## Trajectory Inference
#' ### slingshot
#' Condition Heatmap Cutoffs: p_val_adj <= 0.05 & |log2FC| >= 1  
#' DEGs Heatmap Cutoffs: p_val_adj < 0.05 and if there were > 50 genes, only choose 50 to plot
#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10, results='asis', out.width = "100%"
if (file.exists(paste0(res.loc,"plots/ti"))){
  if(file.exists(paste0(res.loc, "plots/ti/ti_start_smooth.pdf"))){
    
    ReadImageAndTrim(paste0(res.loc, "plots/ti/ti_start_smooth.pdf"))
    
    if(file.exists(paste0(res.loc, "plots/ti/ti_start_smooth.pdf"))){
      ReadImageAndTrim(paste0(res.loc, "plots/ti/ti_deg_between_group.pdf"))
    }
    
    plots <- list.files(paste0(res.loc, "plots/ti/"), full.names=TRUE)
    plots <- plots[grep("slingPseudotime", plots)]
    
    for(i in plots){
      ReadImageAndTrim(i)
    }
  } else if (file.exists(paste0(res.loc,"plots/ti/ti_no_start_not_smooth.pdf"))) {
    ReadImageAndTrim(paste0(res.loc,"plots/ti/ti_no_start_not_smooth.pdf"))
  } else {
    print("slingshot analysis was set to false, if you wanted this analysis please set '-run_sling true'")
  }
} else {
  print("slingshot analysis was set to false, if you wanted this analysis please set '-run_sling true'")
}

#' ### tempora
#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10, results='asis', out.width = "100%"
if(file.exists(paste0(res.loc, "data/se_integrated_tempora_seurat_v3.rds"))){
  plots <- list.files(paste0(res.loc, "plots/ti/"), full.names=TRUE)
  plots <- plots[grep("tempora", plots)]
  
  for(i in plots){
    ReadImageAndTrim(i)
  }
  
  # ReadImageAndTrim(paste0(res.loc, "plots/ti/tempora_screeplot.pdf"))
  # ReadImageAndTrim(paste0(res.loc, "plots/ti/tempora_inferred_lineages.pdf"))
} else {
  print("no time was specified under co-conditions or as the main condition, please rerun with '--co_conditions ...,time,...' or '--main_time true' and ensure your sample names contain it")
}

#' ### psupertime
#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10, results='asis'
if(dir.exists(paste0(res.loc, "plots/ti/psupertime_plots/"))){
  
  plots <- list.files(paste0(res.loc, "plots/ti/psupertime_plots/"), full.names=TRUE, recursive = TRUE)
  plots <- plots[grep("compare_dist", plots)]
  
  for(i in plots){
    ReadImageAndTrim(i)
  }
  
  plots <- list.files(paste0(res.loc, "plots/ti/psupertime_plots/"), full.names=TRUE, recursive = TRUE)
  plots <- plots[grep("psuper_top_20_gene", plots)]
  
  for(i in 1:length(plots)){
    if (i <= 4){
      ReadImageAndTrim(plots[i])
    } else {
      break
    }
  }
} else {
  print("no time was specified under co-conditions or as the main condition, please rerun with '--co_conditions ...,time,...' or '--main_time true' and ensure your sample names contain it")
}

#' ## Differential Abundance analysis with miloR
#' Sorted by character string and numbers, where the first element will have a negative logFC and vice versa. Ex. D100, D56, D70 = D56, D70, D100.  
#' Individual GO pathways in a geyser plot can be found at 
{{if(length(list.files(paste0(res.loc, "plots/da/"), full.names=TRUE)) == 3){"No Differential Abundant neighborhoods were found, so no DEG heatmaps will be made."}}}  
#' Neighborhood cutoff: SpatialFDR < 0.05  
#' DEGs cutoff: |logFC| >= 1 and adj.Pval <= 0.01  
#' GSEA cutoff: p_val_adj <= 0.05
#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10, results='asis'
plots <- list.files(paste0(res.loc, "plots/da/"), full.names=TRUE)
plots <- plots[grep("pval|volcano", plots)]

for(i in plots){
  ReadImageAndTrim(i)
}

plots <- list.files(paste0(res.loc, "plots/da/"), full.names=TRUE)
plots <- plots[grep("umap|fc_dist", plots)]

for(i in plots){
  ReadImageAndTrim(i)
}

#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10, out.width = "50%", fig.show = 'hold', results='asis'
plots <- list.files(paste0(res.loc, "plots/da/"), full.names=TRUE)
plots <- plots[grep("DE_heatmap", plots)]

for(i in plots){
  ReadImageAndTrim(i)
}

plots <- list.files(paste0(res.loc, "plots/da/"), full.names=TRUE)
plots <- plots[grep("gsea", plots)]

for(i in plots){
	  ReadImageAndTrim(i)
}

#' ## Ligand-receptor analysis with CellChat
#' - Bubbleplot is filtered by the missing presence of the communication prob in the other condition, or if there is an abs(log2FC) >= 0.6 -> ~50% increase in prob,  
#' sometimes the dataset cannot be filtered (NA source targets after filtering), so print the unfiltered option in these situations.  
#' - All communication probability graphs can be found in the cellchat_plots folder  
#' - Significant information flow pathways were also plotted separately in the signaling_pathways folder
#+ warning=FALSE, echo=FALSE, fig.height = 10, fig.width = 10, results='asis'
plots <- list.files(paste0(res.loc, "plots/cellchat_plots/"), full.names=TRUE)
#plots <- list.files(plots[1], full.names=TRUE)

ReadImageAndTrim(paste0(plots[1], "/cellchat_interaction_summary_bar.pdf"))

ReadImageAndTrim(paste0(plots[1], "/cellchat_differential_interaction_circle.pdf"))
ReadImageAndTrim(paste0(plots[1], "/cellchat_differential_interaction_heatmap.pdf"))

ReadImageAndTrim(paste0(plots[1], "/cellchat_num_interactions_circle.pdf"))
ReadImageAndTrim(paste0(plots[1], "/cellchat_population_send_receive.pdf"))

ReadImageAndTrim(paste0(plots[1], "/cellchat_compare_outgoing_signal_heatmap.pdf"))
ReadImageAndTrim(paste0(plots[1], "/cellchat_compare_incoming_signal_heatmap.pdf"))
ReadImageAndTrim(paste0(plots[1], "/cellchat_compare_all_signal_heatmap.pdf"))


commun_prob_plots <- list.files(plots[10], full.names=TRUE)

for(i in 1:length(commun_prob_plots)){
  if (i <= 3){
    ReadImageAndTrim(commun_prob_plots[i])
  } else {
    break
  }
}

ReadImageAndTrim(paste0(plots[1], "cellchat_information_flow_compare.pdf"))

if (length(plots) >= 11){
  sig_path_plots <- list.files(plots[11], full.names=TRUE)
  
  if (length(sig_path_plots) == 0){
    print("no significant signalling paths were found")
    
  } else if (length(sig_path_plots) <= 3){
    sig_path_plots_all <- list.files(sig_path_plots, full.names=TRUE, recursive = T)
    for(i in sig_path_plots_all){
      ReadImageAndTrim(i)
    }
    
  } else {
    sig_path_plots_all <- list.files(sig_path_plots[1:3], full.names=TRUE, recursive = T)
    for(i in sig_path_plots_all){
      ReadImageAndTrim(i)
    }
  }
}