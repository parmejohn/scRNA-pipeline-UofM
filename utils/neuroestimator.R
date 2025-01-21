library(neuroestimator)

NeuroestimatorResults <- function(counts_matrix, species){ 
  
  Sys.setenv(RETICULATE_MINICONDA_PATH="/opt/conda")
  reticulate::use_condaenv("neuroestimator")
  
  counts_matrix <- readRDS(counts_matrix)
  #counts_matrix <- as.matrix(se.integrated[["RNA"]]@counts)
  res <- neuroestimator(counts_matrix, species=species) # rmbr to change species if needed
  write.table(res, "neuroestimator_results.txt", quote = FALSE, row.names = T, sep = "\t", col.names = T)
}
