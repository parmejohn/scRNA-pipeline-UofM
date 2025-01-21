#!/opt/R/4.3.2/bin/Rscript
    
library(Seurat)
library(SeuratObject)

set.seed(333)
## set command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

#stop the script if no command line argument
if(length(args)==0){
  print("Please include neuroestimator results")
  stop("Requires command line argument.")
}
  

# args[1] = se_integrated dataset
se.integrated <- readRDS(args[1])
se.integrated[["RNA3"]] <- as(object = se.integrated[["RNA"]], Class = "Assay")
se.integrated <- RenameAssays(se.integrated, "RNA", "RNA5")
se.integrated <- RenameAssays(se.integrated, "RNA3", "RNA")
DefaultAssay(object = se.integrated) <- "RNA"
se.integrated[["RNA5"]] <- NULL
    
saveRDS(se.integrated, "se_integrated_v3.rds")

counts_matrix <- as.matrix(se.integrated[["RNA"]]@counts)
saveRDS(counts_matrix, "se_integrated_v3_counts_matrix.rds")

