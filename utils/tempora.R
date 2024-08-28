TemporaHelper <- function(se.integrated.tempora, gmt_file, group_name = ""){
  
  se.integrated.tempora@cluster.metadata[["label"]] <- gsub( "(.*)-(.*)", "\\2", se.integrated.tempora@cluster.metadata[["label"]]) # remove the cluster prefix
  # get the pathway enrichment profiles from GSVA
  # PCA is used on the different paths that the clusters get in order to remove redundant pathways
  se.integrated.tempora <- CalculatePWProfiles_updated(se.integrated.tempora, 
                                                       gmt_path = gmt_file,
                                                       method="gsva", min.sz = 5, max.sz = 200, parallel.sz = 1)
  pdf(paste0('tempora_screeplot', group_name ,'.pdf'), width = 8, height = 6)
  screeplot(se.integrated.tempora@cluster.pathways.dr, npcs = 25, type = "lines", 
            main = paste0("PCA on pathway enrichment analysis result", group_name))
  graphics.off()
  
  # calculate the scree plot and then find the points with the most differences
  var_explained = se.integrated.tempora@cluster.pathways.dr[["sdev"]]^2 / sum(se.integrated.tempora@cluster.pathways.dr[["sdev"]]^2)
  opt.pc <- which(min(diff(var_explained)) == diff(var_explained)) + 1
  
  # builids trajectory based on clusters pathway enrichment profiles
  se.integrated.tempora <- BuildTrajectory_fixed(se.integrated.tempora, n_pcs = opt.pc, difference_threshold = 0.01)
  
  pdf(paste0('tempora_inferred_lineages', group_name, '.pdf'), width = 8, height = 6)
  se.integrated.tempora <- PlotTrajectory_font_fix(se.integrated.tempora)
  # p2 <- recordPlot()
  # PrintSave(p2, "tempora_inferred_lineages.pdf")
  title(paste0("Time-Series Trajectory Inference", group_name))
  graphics.off()
  
  return(se.integrated.tempora)
}

RunTempora <- function(se.integrated, main_time, species){

  se.integrated[["RNA3"]] <- as(object = se.integrated[["RNA"]], Class = "Assay") #Tempora only works with older seurat objects
  se.integrated$tempora.labels <- Idents(se.integrated)
  
  gmt_file <- DownloadGMT(species)
  
  cluster.names <- NA
  se.integrated@assays[["RNA3_var"]] <- CreateAssayObject(GetAssayData(se.integrated, assay = "RNA3")[se.integrated@assays[["RNA3"]]@var.features,]) # take the 2000 most variable features
  if (all(levels(se.integrated@active.ident) != levels(se.integrated$seurat_clusters))){
    if (main_time == "no"){
      if (length(unique(se.integrated@meta.data[["group"]])) >= 2){
        for (i in unique(se.integrated@meta.data[["group"]])){
          se.integrated.filt <- subset(se.integrated, group == i)
          se.integrated.tempora <-  ImportSeuratObject(se.integrated.filt,  clusters = "tempora.labels", 
                                                       timepoints = "time", assayType = "RNA3_var", 
                                                       assaySlot = "data", 
                                                       cluster_labels = levels(se.integrated$tempora.labels),
                                                       timepoint_order = mixedsort(unique(se.integrated@meta.data[["time"]])))
          se.integrated.tempora <- TemporaHelper(se.integrated.tempora, gmt_file, paste0("_", i))
        }
      } else {
        se.integrated.tempora <-  ImportSeuratObject(se.integrated,  clusters = "tempora.labels", 
                                                     timepoints = "time", assayType = "RNA3_var", 
                                                     assaySlot = "data", 
                                                     cluster_labels = levels(se.integrated$tempora.labels),
                                                     timepoint_order = mixedsort(unique(se.integrated@meta.data[["time"]])))
        se.integrated.tempora <- TemporaHelper(se.integrated.tempora, gmt_file)
      }
    } else if (main_time == "yes") {
      se.integrated.tempora <-  ImportSeuratObject(se.integrated,  clusters = "tempora.labels", 
                                                   timepoints = "group", assayType = "RNA3_var", 
                                                   assaySlot = "data", 
                                                   cluster_labels = levels(se.integrated$tempora.labels),
                                                   timepoint_order = mixedsort(unique(se.integrated@meta.data[["group"]])))
      se.integrated.tempora <- TemporaHelper(se.integrated.tempora, gmt_file)
    }
  } else {
    se.integrated.tempora <-  ImportSeuratObject(se.integrated,  clusters = "seurat_clusters", 
                                                 timepoints = "time", assayType = "RNA3_var", 
                                                 assaySlot = "data", 
                                                 timepoint_order = mixedsort(unique(se.integrated@meta.data[["time"]])))
    se.integrated.tempora <- TemporaHelper(se.integrated.tempora, gmt_file)
  }
  return(se.integrated.tempora)
}


