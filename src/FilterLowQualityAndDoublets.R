#!/usr/local/bin/Rscript

library(argparse)
library(Seurat)
library(ggplot2)
library(scDblFinder)
library(SingleCellExperiment)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(tidyverse)
library(Signac)

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
parser$add_argument(
  '-species',
  '--s',
  type = "character",
  required = TRUE,
  nargs = 1,
  help = 'Species name (Mus musculus, Homo sapiens); CASE-SENSITIVE'
)
parser$add_argument(
  '-coconditions',
  type = "character",
  required = TRUE,
  nargs = '*',
  help = 'Co-conditions listed'
)
parser$add_argument(
  '-atac',
  type = "character",
  required = TRUE,
  nargs = 1,
  help = 'Is this a multiome experiment?'
)
parser$add_argument(
  '-original_files',
  type = "character",
  required = TRUE,
  nargs = 1,
  help = 'Contains CellRanger count outputs folder seperated by condition'
)
args <- parser$parse_args()

indir <- args$i

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


##### Loading data as Seurat Objects #####
print("Loading data as Seurat Objects")
foldernames <-
  list.dirs(path = indir,
            full.names = TRUE,
            recursive = T)
foldernames <- foldernames[grepl("soupx", foldernames)]

print(indir)
se.list <-
  lapply(foldernames, Read10X) %>% lapply(CreateSeuratObject) # load the SoupX corrected data as a adj.matrix and then convert it into a SeuratObject
print(se.list)

species <- args$s

if (args$atac == "yes"){
  if (species == "musmusculus"){
    
    
    annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
    
  } else if (species == "homosapiens"){
    
    
    annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
    
  }
  # find common peak set
  bedfiles <-
    list.files(
      path = args$original_files,
      pattern = "atac_peaks.bed",
      full.names = TRUE,
      recursive = T,
      include.dirs = T
    )
  print("read tables")
  bed.list <- lapply(bedfiles, read.table, header = FALSE, col.names = c("chr", "start", "end"), sep = "\t")
  print("make granges")
  gr.list <- lapply(bed.list, makeGRangesFromDataFrame)
  #gr.list <- GRangesList(gr.list)
  gr.vector <- unlist(gr.list)
  print(class(gr.vector))
  print("combining peaks")
  #combined.peaks <- Signac::reduce(x = gr.vector)
  combined.peaks <- UnifyPeaks(object.list = gr.list)
  peakwidths <- width(combined.peaks)
  combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20] # filtering bad peaks based on length
  
  # should be the same length and order as the se.list below
  filenames <-
    list.files(
      path = args$original_files,
      pattern = "filtered_feature_bc_matrix.h5|fragments.tsv.gz$",
      full.names = TRUE,
      recursive = T,
      include.dirs = T
    )
  print(filenames)
  filenames <- filenames[grepl("outs", filenames)]
  
  list_of_pairs <- list()
  for (i in 1:length(filenames)) {
    if (i %% 2 == 1) {
      temp_list <- list(c(filenames[i], filenames[i + 1]))
      list_of_pairs <- append(list_of_pairs, temp_list)
    }
  }
  
  for (i in 1:length(se.list)){
    filt.matrix <- Read10X_h5(list_of_pairs[[i]][2], use.names = T)
    
    se.temp <- se.list[[i]]
    
    existing_cells <- colnames(se.temp)
    #filtered_chromatin_data <- filt.matrix$Peaks[, existing_cells, drop = FALSE]
    
    frags.obj <- CreateFragmentObject(path = list_of_pairs[[i]][1])
    
    feature.counts <- FeatureMatrix(
      fragments = frags.obj,
      features = combined.peaks
    )
    
    filtered.feature.counts <- feature.counts[, existing_cells, drop = FALSE]
    
    ## need to reload the chromatin assay after
    chrom_assay <- CreateChromatinAssay(
      counts = filtered.feature.counts,
      sep = c(":", "-"),
      fragments = frags.obj,
      annotation = annotation
    )
    
    se.temp[["ATAC"]] <- chrom_assay
    DefaultAssay(se.temp) <- "RNA"
    se.list[[i]] <-  se.temp
  }
}

# set up sample and conditions
print("Setting up se objects")
for (i in 1:length(se.list)) {
  se.list[[i]]@misc <-
    list(sub(".*\\/(.*)", "\\1", foldernames[i])) # save sample folder name in SeuratObject metadata
  se.list[[i]] <-
    SetIdent(se.list[[i]], value = sub(".*\\/(.*)", "\\1", foldernames[i]))
  sample.name <- sub(".*\\/(.*)", "\\1", foldernames[i])
  if (!grepl("_", sample.name)){
    stop("there should be at least 1 underscore in your sample names, please relabel them")
  }
  se.list[[i]]$sample <- sample.name
  se.list[[i]]$group <-
    sub(".*\\/(.*)\\/.*", "\\1", foldernames[i]) # save condition folder name under group
  if (args$coconditions[1] != 'none'){
    if (str_count(sample.name, "_") < length(args$coconditions) + 2){
      stop("you might have forgetten an underscore somewhere when trying to seperate your conditions, please check again")
    }
    for (j in 1:length(args$coconditions)){
      se.list[[i]] <- AddMetaData(se.list[[i]], sapply(strsplit(sample.name, "_"), function(x) x[j+2]), args$coconditions[j])
    }
    Misc(se.list[[i]], slot = "co.conditions") <- args$coconditions
  }
}
saveRDS(se.list, "se_list_raw.rds")

##### Basic QC #####
print("Basic QC")
se.filtered.list <-
  lapply(se.list, BasicQC, species = species, atac = args$atac) # removal of low quality cells by percentage of mitochondrial reads, number of genes, and number of UMIs
saveRDS(se.filtered.list, "se_filtered_list.rds")

##### Doublet Removal #####
print("Doublet removal")
se.filtered.doublets.list <- lapply(se.filtered.list, DoubletQC, atac = args$atac)
saveRDS(se.filtered.doublets.list, "se_filtered_doublets_list.rds")

se.filtered.singlets.list <- lapply(se.filtered.doublets.list, subset, subset = scDblFinder.class  == "singlet")

saveRDS(se.filtered.singlets.list, "se_filtered_singlets_list.rds")
