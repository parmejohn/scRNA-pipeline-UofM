RunTempora <- function(se.integrated){

  se.integrated[["RNA3"]] <- as(object = se.integrated[["RNA"]], Class = "Assay") #Tempora only works with older seurat objects
  se.integrated$tempora.labels <- Idents(se.integrated)
  
  cluster.names <- NA
  if (all(levels(se.integrated@active.ident) != levels(se.integrated$seurat_clusters))){
    
    ## This is to grab which cluster name correlates with which seurat clustering labels
    # prediction.scores <- se.integrated@meta.data[, grepl("^prediction|RNA_snn_res.1", names(se.integrated@meta.data))]
    # prediction.scores <- prediction.scores[,-which(names(prediction.scores) == "prediction.score.max")]
    # 
    # colnames(prediction.scores) <- gsub("prediction.score.", "", colnames(prediction.scores))
    # prediction.scores <- melt(prediction.scores, id.vars = "RNA_snn_res.1", variable.name = "source", value.name = "score")
    # 
    # prediction.matrix <- tapply(prediction.scores$score, list(prediction.scores$RNA_snn_res.1, prediction.scores$source), median)
    # cluster.names <- colnames(prediction.matrix)[max.col(prediction.matrix,ties.method="first")]
    # 
    # se.integrated <- SetIdent(se.integrated, value = se.integrated@meta.data$seurat_clusters)
    se.integrated@assays[["RNA3_var"]] <- CreateAssayObject(GetAssayData(se.integrated, assay = "RNA3")[se.integrated@assays[["RNA3"]]@var.features,]) # take the 2000 most variable features
    # 
    # se.integrated.tempora <-  ImportSeuratObject(se.integrated,  clusters = "seurat_clusters", 
    #                                              timepoints = "time", assayType = "RNA3_var", 
    #                                              assaySlot = "data", 
    #                                              cluster_labels = cluster.names,
    #                                              timepoint_order = mixedsort(unique(se.integrated@meta.data[["time"]])))
    
    se.integrated.tempora <-  ImportSeuratObject(se.integrated,  clusters = "tempora.labels", 
                                                 timepoints = "time", assayType = "RNA3_var", 
                                                 assaySlot = "data", 
                                                 cluster_labels = levels(se.integrated$tempora.labels),
                                                 timepoint_order = mixedsort(unique(se.integrated@meta.data[["time"]])))
  } else {
    se.integrated.tempora <-  ImportSeuratObject(se.integrated,  clusters = "seurat_clusters", 
                                                 timepoints = "time", assayType = "RNA3_var", 
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
  se.integrated.tempora <- BuildTrajectory_fixed(se.integrated.tempora, n_pcs = opt.pc, difference_threshold = 0.01)
  
  # font_import(prompt = F)
  # #font_import(pattern="Arial", prompt = F)
  # fonts()
  # loadfonts()

  pdf('tempora_inferred_lineages.pdf', width = 8, height = 6)
  se.integrated.tempora <- PlotTrajectory_font_fix(se.integrated.tempora, hgap=1, vgap=1)
  # p2 <- recordPlot()
  # PrintSave(p2, "tempora_inferred_lineages.pdf")
  title("Time-Series Trajectory Inference")
  graphics.off()
  
  return(se.integrated.tempora)
}


