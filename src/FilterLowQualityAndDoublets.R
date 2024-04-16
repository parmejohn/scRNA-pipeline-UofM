#!/usr/local/bin/Rscript
#source(paste0(dirname(dirname(dirname(getwd()))),"/utils/qc.R"))

library(argparse)
library(Seurat)
library(ggplot2)
library(scDblFinder)
library(SingleCellExperiment)
library(dplyr)
library(tidyverse)

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
#filenames.filt <- filenames[grepl("filtered", filenames)]

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
  }
  Misc(se.list[[i]], slot = "co.conditions") <- args$coconditions
}
saveRDS(se.list, "se_list_raw.rds")

##### Basic QC #####
print("Basic QC")
se.filtered.list <-
  lapply(se.list, BasicQC, species = args$species) # removal of low quality cells by percentage of mitochondrial reads, number of genes, and number of UMIs
saveRDS(se.filtered.list, "se_filtered_list.rds")

##### Doublet Removal #####
print("Doublet removal")
se.filtered.doublets.list <- lapply(se.filtered.list, DoubletQC)
saveRDS(se.filtered.doublets.list, "se_filtered_doublets_list.rds")

se.filtered.singlets.list <- lapply(se.filtered.doublets.list, subset, subset = scDblFinder.class  == "singlet")

saveRDS(se.filtered.singlets.list, "se_filtered_singlets_list.rds")
