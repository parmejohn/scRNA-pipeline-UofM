#Estimate pathway enrichment profiles of clusters
CalculatePWProfiles_updated <- function (object, gmt_path, method = "gsva", min.sz = 5, max.sz = 200, 
                                         parallel.sz = 1) 
{
  if (class(object)[1] != "Tempora") {
    stop("Not a valid Tempora object")
  }
  cat("Calculating cluster average gene expression profile...")
  #exprMatrix <- object@data
  exprMatrix <- Matrix(object@data, sparse = TRUE)
  exprMatrix_bycluster <- list()
  pathwaygmt <- GSEABase::getGmt(gmt_path)
  for (i in sort(unique(object@meta.data$Clusters))) {
    exprMatrix_bycluster[[i]] <- rowMeans(exprMatrix[, which(colnames(exprMatrix) %in% 
                                                               rownames(object@meta.data)[which(object@meta.data$Clusters == 
                                                                                                  i)])])
  }
  names(exprMatrix_bycluster) <- sort(unique(object@meta.data$Clusters))
  exprMatrix_bycluster <- Matrix(do.call(cbind, exprMatrix_bycluster), sparse = TRUE)
  colnames(exprMatrix_bycluster) <- sort(unique(object@meta.data$Clusters))
  rownames(exprMatrix_bycluster) <- rownames(exprMatrix)
  # GSVA was updated to require this gsvaParam arguement instead of the default method that Tempora was using
  print("Calculating GSVA Param")
  gsvapar <- GSVA::gsvaParam(exprMatrix_bycluster, pathwaygmt, minSize = min.sz, maxSize = max.sz)
  cat("\nCalculating cluster pathway enrichment profiles...\n")
  gsva_bycluster <- GSVA::gsva(gsvapar, BPPARAM = SerialParam(progressbar = TRUE))
  colnames(gsva_bycluster) <- colnames(exprMatrix_bycluster)
  object@cluster.pathways <- gsva_bycluster
  gsva_bycluster_pca <- prcomp(t(gsva_bycluster), scale = T, 
                               center = T)
  #screeplot(gsva_bycluster_pca, npcs = 25, type = "lines", 
  #          main = "PCA on pathway enrichment analysis result")
  object@cluster.pathways.dr <- gsva_bycluster_pca
  validObject(object)
  return(object)
}

PlotTrajectory_font_fix <- function (object, layout = NULL, ...) 
{
  if (class(object)[1] != "Tempora") {
    stop("Not a valid Tempora object")
  }
  if (is.null(object@trajectory)) {
    stop("BuildTrajectory has not been run. See ?Tempora::BuildTrajectory for details")
  }
  edge_graph <- igraph::graph_from_data_frame(d = object@trajectory, 
                                              vertices = object@cluster.metadata, directed = T)
  if (is.null(layout)) {
    l <- igraph::layout_with_sugiyama(edge_graph, layers = object@cluster.metadata$Cluster_time_score, 
                                      maxiter = 1000)
    if (length(levels(object@meta.data$Timepoints)) > 9) {
      colours <- colorRampPalette(RColorBrewer::brewer.pal(7, 
                                                           "YlOrRd"))
      plot.igraph(edge_graph, ylim = c(-1, 1), layout = l$layout, 
                  ylab = "Inferred time", vertex.shape = "pie", 
                  vertex.pie = lapply(1:nrow(object@cluster.metadata), 
                                      function(x) as.numeric(object@cluster.metadata[x, 
                                                                                     2:((length(levels(object@meta.data$Timepoints))) + 
                                                                                          1)])), vertex.pie.color = list(colours(length(levels(object@meta.data$Timepoints)))), 
                  pie.border = list(rep("white", 4)), vertex.frame.color = "white", 
                  edge.arrow.size = 0.5, edge.width = 1.5,
                  vertex.label.color = "black", edge.lty = E(edge_graph)$type, 
                  ...)
      axis(side = 2, at = c(-1, 1), labels = c("Late", 
                                               "Early"), las = 1)
      legend("topleft", legend = levels(object@meta.data$Timepoints), 
             fill = colours, bty = "n", border = "black")
    }
    else {
      colours <- brewer.pal(length(levels(object@meta.data$Timepoints)), 
                            "YlOrRd")
      plot.igraph(edge_graph, ylim = c(-1, 1), ylab = "Inferred time", 
                  layout = l$layout, vertex.shape = "pie", vertex.pie = lapply(1:nrow(object@cluster.metadata), 
                                                                               function(x) as.numeric(object@cluster.metadata[x, 
                                                                                                                              2:((length(levels(object@meta.data$Timepoints))) + 
                                                                                                                                   1)])), vertex.pie.color = list(colours), 
                  pie.border = list(rep("white", length(levels(object@meta.data$Timepoints)))), 
                  vertex.frame.color = "white",
                  vertex.label.color = "black", edge.lty = E(edge_graph)$type, 
                  ...)
      legend("topleft", legend = levels(object@meta.data$Timepoints), 
             fill = colours, bty = "n", border = "black")
      axis(side = 2, at = c(-1, 1), labels = c("Late", 
                                               "Early"), las = 1)
    }
    object@layouts <- l$layout
  }
  else {
    if (length(levels(object@meta.data$Timepoints)) > 9) {
      colours <- colorRampPalette(RColorBrewer::brewer.pal(7, 
                                                           "YlOrRd"))
      plot.igraph(edge_graph, ylim = c(-1, 1), layout = layout, 
                  ylab = "Inferred time", vertex.shape = "pie", 
                  vertex.pie = lapply(1:nrow(object@cluster.metadata), 
                                      function(x) as.numeric(object@cluster.metadata[x, 
                                                                                     2:((length(levels(object@meta.data$Timepoints))) + 
                                                                                          1)])), vertex.pie.color = list(colours(length(levels(object@meta.data$Timepoints)))), 
                  pie.border = list(rep("white", 4)), vertex.frame.color = "white", 
                  edge.arrow.size = 0.5, edge.width = 1.5,
                  vertex.label.color = "black", edge.lty = E(edge_graph)$type, 
                  ...)
      axis(side = 2, at = c(-1, 1), labels = c("Late", 
                                               "Early"), las = 1)
      legend("topleft", legend = levels(object@meta.data$Timepoints), 
             fill = colours, bty = "n", border = "black")
    }
    else {
      colours <- brewer.pal(length(levels(object@meta.data$Timepoints)), 
                            "YlOrRd")
      plot.igraph(edge_graph, ylim = c(-1, 1), ylab = "Inferred time", 
                  layout = layout, vertex.shape = "pie", vertex.pie = lapply(1:nrow(object@cluster.metadata), 
                                                                             function(x) as.numeric(object@cluster.metadata[x, 
                                                                                                                            2:((length(levels(object@meta.data$Timepoints))) + 
                                                                                                                                 1)])), vertex.pie.color = list(colours), 
                  pie.border = list(rep("white", length(levels(object@meta.data$Timepoints)))), 
                  vertex.frame.color = "white",
                  vertex.label.color = "black", edge.lty = E(edge_graph)$type, 
                  ...)
      legend("topleft", legend = levels(object@meta.data$Timepoints), 
             fill = colours, bty = "n", border = "black")
      axis(side = 2, at = c(-1, 1), labels = c("Late", 
                                               "Early"), las = 1)
    }
  }
  validObject(object)
  return(object)
}


DownloadGMT <- function(){
  options(timeout=600)
  
  ##### Downloading Gene Set File #####
  gmt_url = "http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/"
  
  #list all the files on the server
  filenames = getURL(gmt_url)
  tc = textConnection(filenames)
  contents = readLines(tc)
  close(tc)
  
  #get the gmt that has all the pathways and does not include terms inferred from electronic annotations(IEA)
  #start with gmt file that has pathways only
  
  rx = gregexpr("(?<=<a href=\")(.*.GOBP_AllPathways_noPFOCR_no_GO_iea.*.)(.gmt)(?=\">)",
                contents, perl = TRUE)
  
  gmt_file = unlist(regmatches(contents, rx))
  dest_gmt_file <- file.path(getwd(),gmt_file )
  download.file(
    paste(gmt_url,gmt_file,sep=""),
    destfile=dest_gmt_file
  )
  return(gmt_file)
}

RunTempora <- function(se.integrated){

  se.integrated <- readRDS(input)
  se.integrated[["RNA3"]] <- as(object = se.integrated[["RNA"]], Class = "Assay") #Tempora only works with older seurat objects
  
  cluster.names <- NA
  if (all(levels(se.integrated@active.ident) != levels(se.integrated$seurat_clusters))){
    
    ## This is to grab which cluster name correlates with which seurat clustering labels
    prediction.scores <- se.integrated@meta.data[, grepl("^prediction|RNA_snn_res.1", names(se.integrated@meta.data))]
    prediction.scores <- prediction.scores[,-which(names(prediction.scores) == "prediction.score.max")]
    
    colnames(prediction.scores) <- gsub("prediction.score.", "", colnames(prediction.scores))
    prediction.scores <- melt(prediction.scores, id.vars = "RNA_snn_res.1", variable.name = "source", value.name = "score")
    
    prediction.matrix <- tapply(prediction.scores$score, list(prediction.scores$RNA_snn_res.1, prediction.scores$source), median)
    cluster.names <- colnames(prediction.matrix)[max.col(prediction.matrix,ties.method="first")]
    
    se.integrated <- SetIdent(se.integrated, value = se.integrated@meta.data$seurat_clusters)
    se.integrated.tempora <-  ImportSeuratObject(se.integrated,  clusters = "seurat_clusters", 
                                                 timepoints = "time", assayType = "RNA3", 
                                                 assaySlot = "data", 
                                                 cluster_labels = cluster.names,
                                                 timepoint_order = mixedsort(unique(se.integrated@meta.data[["time"]])))
  } else {
    se.integrated.tempora <-  ImportSeuratObject(se.integrated,  clusters = "seurat_clusters", 
                                                 timepoints = "time", assayType = "RNA3", 
                                                 assaySlot = "data", 
                                                 timepoint_order = mixedsort(unique(se.integrated@meta.data[["time"]])))
  }
  
  gmt_file <- DownloadGMT()
  
  # get the pathway enrichment profiles from GSVA
  # PCA is used on the different paths that the clusters get in order to remove redundant pathways
  se.integrated.tempora <- CalculatePWProfiles_updated(se.integrated.tempora, 
                                                       gmt_path = gmt_file,
                                                       method="gsva", min.sz = 5, max.sz = 200, parallel.sz = 1)
  pdf('tempora_screeplot.pdf', width = 8, height = 6)
  screeplot(se.integrated.tempora@cluster.pathways.dr, npcs = 25, type = "lines", 
            main = "PCA on pathway enrichment analysis result")
  # p1 <- recordPlot()
  # PrintSave(p1, "tempora_screeplot.pdf")
  graphics.off()
  
  # calculate the scree plot and then find the points with the most differences
  var_explained = se.integrated.tempora@cluster.pathways.dr[["sdev"]]^2 / sum(se.integrated.tempora@cluster.pathways.dr[["sdev"]]^2)
  opt.pc <- which(min(diff(var_explained)) == diff(var_explained)) + 1
  
  # builids trajectory based on clusters pathway enrichment profiles
  se.integrated.tempora <- BuildTrajectory(se.integrated.tempora, n_pcs = opt.pc, difference_threshold = 0.01)
  
  # font_import(prompt = F)
  # #font_import(pattern="Arial", prompt = F)
  # fonts()
  # loadfonts()

  pdf('tempora_inferred_lineages.pdf', width = 8, height = 6)
  se.integrated.tempora <- PlotTrajectory_font_fix(se.integrated.tempora)
  # p2 <- recordPlot()
  # PrintSave(p2, "tempora_inferred_lineages.pdf")
  title("Time-Series Trajectory Inference")
  graphics.off()
  
  return(se.integrated.tempora)
}