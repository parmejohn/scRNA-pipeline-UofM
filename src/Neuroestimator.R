#!/opt/conda/envs/neuroestimator/bin/Rscript
library(neuroestimator)

set.seed(333)

## set command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

#stop the script if no command line argument
if(length(args)==0){
  print("Please include seurat object")
  stop("Requires command line argument.")
}


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
))), "/utils/neuroestimator.R"))

if (args[2] == 'homosapiens') {
  species <- 'hsapiens'
} else if (args[2] == 'musmusculus'){
  species <- 'mmusculus'
}

# args[1] = gene expression table path
# args[2] = species
NeuroestimatorResults(args[1], species)
