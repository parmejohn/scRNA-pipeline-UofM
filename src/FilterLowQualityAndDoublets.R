#!/usr/local/bin/Rscript
source(paste0(dirname(dirname(dirname(getwd()))),"/utils/qc.R"))

library(argparse)
library(Seurat)
library(ggplot2)
library(scDblFinder)
library(SingleCellExperiment)
library(dplyr)
library(tidyverse)

set.seed(333)

parser <- ArgumentParser(description='Process scRNA-seq data, while performing comparitive analyses')
parser$add_argument('-indir', '--i',  type="character", required=TRUE, nargs=1, help='Contains CellRanger count outputs folder seperated by condition')
args <- parser$parse_args()

indir <- args$i

##### Loading data as Seurat Objects #####
print("Loading data as Seurat Objects")
foldernames <- list.dirs(path = paste0(args$i), full.names = TRUE, recursive = T)
foldernames <- foldernames[grepl("soupx", foldernames)]

se.list <- lapply(foldernames, Read10X) %>% lapply(CreateSeuratObject) # load the SoupX corrected data as a adj.matrix and then convert it into a SeuratObject

#filenames.filt <- filenames[grepl("filtered", filenames)]

# set up sample and conditions
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
