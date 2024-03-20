#!/usr/bin/env Rscript
set.seed(333) # ensure that results are reproducible when using the same data

#### Libraries ####
#library(argparse)
packages <- c("argparse", "dplyr", "ggplot2", "ggpubr", "ggrepel", "Seurat", "SeuratObject",
              "SingleCellExperiment", "SoupX", "scDblFinder", "tibble", "tidyverse",
              "org.Hs.eg.db", "hdf5r", "DropletUtils", "SeuratWrappers", "presto", "scCustomize",
              "pheatmap", "msigdbr", "fgsea", "slingshot", "tradeSeq", "gam", "ComplexHeatmap",
              "RColorBrewer", "escape", "NbClust", "glmGamPoi")

invisible(lapply(packages, library, character.only = TRUE))

#### Setting up argparse ####
parser <- ArgumentParser(description='Process scRNA-seq data, while performing comparitive analyses')

# required
parser$add_argument('-indir', '--i',  type="character", required=TRUE, nargs=1, help='Contains CellRanger count outputs folder seperated by condition')
parser$add_argument('-outdir', '--o', type="character", required=TRUE, nargs=1, help='Folder to create sc_analysis folder')
parser$add_argument('-species', '--s',type="character", required=TRUE, nargs=1, help='Species name (Mus musculus, Homo sapiens); CASE-SENSITIVE')

# optional
parser$add_argument('-reference_seurat', type="character", nargs=1, help='Reference Seurat(s) object file path with pre-labelled clusters') # need to change to many options later
parser$add_argument('-clusters_optimal', type="integer", nargs=1, help='Optimal clusters for dimensional reductions and clustering algorithms')
parser$add_argument('-resolution', type="integer", nargs=1, help='Change resolution for Seurat UMAP, higher values will lead to more clustering')
parser$add_argument('-beginning_cluster', type="character", nargs=1, help='Beginning cluster for trajectory inference with slingshot')


args <- parser$parse_args()

# print(args)
print(paste0(args$i))
print(paste0(args$o))
print(paste0(args$s))
# print(paste0(args$reference_seurat))
# print(paste0(args$clusters_optimal))

 # if (length(args$clusters_optimal) == 0) {
 #   print('yes')
 # }
 
#### Variables ####
#link to folder which is separated by conditions, has sample count matrix folders inside
#these condition names will be used downstream, so make sure they what is wanted
indir <- args$i

# set output folder; probably dont want in same folder as the indir
outdir <- args$o

# create folder outputs
dir.create(paste(outdir, '/analysis/', sep=''))
dir.create(paste(outdir, '/analysis/data/', sep=''))
dir.create(paste(outdir, '/analysis/data/qc/', sep=''))
dir.create(paste(outdir, '/analysis/plots/', sep=''))

plotpdf <- paste(outdir, '/analysis/plots/', sep='')
outdata <- paste(outdir, '/analysis/data/', sep='')

# parameters
opt.clusters <- 0

#### Functions ####
#current working directory should be the downloaded folder from github
file.sources = list.files("./utils/", pattern="*.R",
                          full.names=TRUE, ignore.case=TRUE)
sapply(file.sources, source, .GlobalEnv)

#### Pre-analysis QC ####
##### SoupX filtering #####
print("SoupX filtering")
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

# perform Ambient RNA correction; see function in qc.R script for more info
# files will be saved under second arg
lapply(list_of_pairs, AmbientRNARemoval, outdir)

##### Loading data as Seurat Objects #####
print("Loading data as Seurat Objects")
foldernames <- list.dirs(path = paste(outdir, "analysis/data/qc/", sep=''), full.names = TRUE, recursive = T)
foldernames <- foldernames[grepl("soupx", foldernames)]

se.list <- lapply(foldernames, Read10X) %>% lapply(CreateSeuratObject) # load the SoupX corrected data as a adj.matrix and then convert it into a SeuratObject

#filenames.filt <- filenames[grepl("filtered", filenames)]

for (i in 1:length(se.list)){
  se.list[[i]]@misc <- list(sub(".*\\/(.*)", "\\1", foldernames[i])) # save sample folder name in SeuratObject metadata
  se.list[[i]] <- SetIdent(se.list[[i]], value = sub(".*\\/(.*)", "\\1", foldernames[i]))
  se.list[[i]]$sample <- sub(".*\\/(.*)", "\\1", foldernames[i])
  se.list[[i]]$group <- sub(".*\\/(.*)\\/.*", "\\1", foldernames[i]) # save condition folder name under group
}

##### Basic QC #####
print("Basic QC")
se.filtered.list <- lapply(se.list, BasicQC, plotpdf) # removal of low quality cells by percentage of mitochondrial reads, number of genes, and number of UMIs

##### Doublet Removal #####
print("Doublet removal")
se.filtered.singlets.list <- lapply(se.filtered.list, DoubletQC)

#### Integrating ####
print("Integrating")
# Using CCA for batch correction; https://hbctraining.github.io/scRNA-seq_online/lessons/06_integration.html
se.integrated <- IntegrateSamples(se.filtered.singlets.list, group)

#### Dimensional Reduction ####
print("Dimensional reduction")
if (length(args$clusters_optimal) == 0) {
  opt.clusters <- NbClust(se.integrated@reductions[["integrated.cca"]]@feature.loadings, distance = "euclidean", min.nc=10, max.nc = 20,
                          method = "complete", index = "ch")$Best.nc[1] %>% unname()
} else {
  opt.clusters <- args$clusters_optimal
}

integrated_elbow <- ElbowPlot(se.integrated) # have to use visual check to find optimal # of clusters for now
# this can always be increased if want to seperate clusters more ie. finding rare cell populations
PrintSave(integrated_elbow, "integrated_elbow_plot.pdf", path=plotpdf)

if (length(args$resolution) == 0) {
  se.integrated <- SeuratDimReduction(se.integrated, 1:opt.clusters, plotpdf, 'group')
} else {
  se.integrated <- SeuratDimReduction(se.integrated, 1:opt.clusters, plotpdf, 'group', args$resolution)
}



#se.integrated <- readRDS("/home/projects/sc_pipelines/analysis/data/se_integrated.rds")

#### Identifying Cell Markers ####
print("IdentifyCellMarkers")
IdentifyCellMarkers(se.integrated, outdata, plotpdf)

##### external dataset #####
print("External mapping")
#se.michalski <- readRDS(paste0(outdata, "/se_michalski.rds"))

if (length(args$reference_seurat) != 0) {
  se.reference <- readRDS(args$reference_seurat)
  se.integrated <- ReferenceMarkerMapping(se.reference, se.integrated, opt.clusters, plotpdf) # active.ident should correlate to the cluster identification
  saveRDS(se.integrated, "/home/projects/sc_pipelines/analysis/data/se_integrated_auto_label.rds")
}


#### Post-identification analysis ####
##### Trajectory inference #####
print("Trajectory Inference")
# https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html#basic_processing_with_seurat_pipeline
if (length(args$beginning_cluster) == 0) {
  TrajectoryInferenceSlingshot(se.integrated, plotpdf, data.path= outdata)
} else {
  TrajectoryInferenceSlingshot(se.integrated, plotpdf, data.path= outdata, start.clus=args$beginning_cluster)
}
  
##### Comparitive DESEQ2 and GSEA #####
print("DESEQ2 and GSEA")
if (length(unique(se.integrated@meta.data[["group"]])) == 2){
  DESeq2ConditionPerCluster(se.integrated, plotpdf, args$s)
}

##### single-cell specific GSEA #####
print("single-cell GSEA")
EscapeGSEA(se.integrated, args$s, plotpdf, outdata)

##### Differential Abundance #####
print("Differential Abundance")
# these packages have conflictions with miloR so they need to be unloaded
unloadNamespace('DropletUtils')
unloadNamespace('scDblFinder')
unloadNamespace('scran')
library(miloR)

# make sure idents are set correctly
Idents(object = se.integrated) <- "celltype"
DifferentialAbundanceMilo(se.integrated, 'sample', 'group', 30, opt.clusters, plotpdf, 'INTEGRATED.CCA')z
