library(neuroestimator)

NeuroestimatorResults <- function(se.integrated, species){ 
  
  Sys.setenv(RETICULATE_MINICONDA_PATH="/opt/conda")
  reticulate::use_condaenv("neuroestimator")
  
  res <- neuroestimator(as.matrix(se.integrated[["RNA"]]@counts), species=species) # rmbr to change species if needed
  write.table(res, "neuroestimator_results.txt", quote = FALSE, row.names = T, sep = "\t", col.names = T)
}
