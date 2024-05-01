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
  object@cluster.metadata$Cluster_time_score_round <- round(object@cluster.metadata$Cluster_time_score)
  if (is.null(layout)) {
    l <- igraph::layout_with_sugiyama(edge_graph, layers = object@cluster.metadata$Cluster_time_score, 
                                      maxiter = 2000 , hgap=1, vgap=1
    ) # need to make layout have additional levels instead of 0 and 1 -> takes the floor of the time score value
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

BuildTrajectory_fixed <- function (object, n_pcs, difference_threshold = 0.01, loadings = 0.4) 
{
  if (class(object)[1] != "Tempora") {
    stop("Not a valid Tempora object")
  }
  if (n_pcs > ncol(object@cluster.pathways.dr$rotation)) {
    stop("Number of PCs selected exceeds number of PCs calculated")
  }
  significant_pathways_list <- gsva_pca <- list()
  for (i in 1:n_pcs) {
    genes_scaled <- scale(object@cluster.pathways.dr$rotation[, 
                                                              i])
    significant_pathways_list[[i]] <- object@cluster.pathways[which(rownames(object@cluster.pathways) %in% 
                                                                      names(which(genes_scaled[, 1] > loadings | genes_scaled[, 
                                                                                                                              1] < (-1 * loadings)))), ]
    gsva_pca[[i]] <- colMeans(significant_pathways_list[[i]])
  }
  gsva_pca <- Reduce(rbind, gsva_pca)
  rownames(gsva_pca) <- paste0("PC", seq(1:nrow(gsva_pca)))
  mi_network <- bnlearn::aracne(as.data.frame(gsva_pca))
  edges_df <- as.data.frame(mi_network$arcs)
  #edges_df$to <- as.numeric(as.character(edges_df$to))
  #edges_df$from <- as.numeric(as.character(edges_df$from))
  edges_df$from_clusterscore <- unlist(sapply(edges_df$from, 
                                              function(x) object@cluster.metadata$Cluster_time_score[object@cluster.metadata$Id == 
                                                                                                       x]))
  edges_df$to_clusterscore <- unlist(sapply(edges_df$to, function(x) object@cluster.metadata$Cluster_time_score[object@cluster.metadata$Id == 
                                                                                                                  x]))
  edges_df$direction <- ifelse((abs(edges_df$to_clusterscore - 
                                      edges_df$from_clusterscore)/(0.5 * (edges_df$to_clusterscore + 
                                                                            edges_df$from_clusterscore))) < difference_threshold, 
                               "bidirectional", "unidirectional")
  edges_df <- edges_df[-which(edges_df$from_clusterscore > 
                                edges_df$to_clusterscore), ]
  # edges_df$id <- ifelse(as.numeric(edges_df$from) > as.numeric(edges_df$to), 
  #                       paste0(edges_df$from, edges_df$to), paste0(edges_df$to, 
  #                                                                  edges_df$from))
  edges_df$id <-paste0(edges_df$to, edges_df$from)
  edges_df <- edges_df[!duplicated(edges_df$id), ]
  edges_df <- edges_df[, -6]
  edges_df$type <- ifelse(edges_df$direction == "bidirectional",
                          3, 1)
  object@trajectory <- edges_df
  object@n.pcs <- n_pcs
  return(object)
}