set.seed(333)

EscapeGSEA <- function(se.integrated, species, pathways, plots.format){
  se.integrated$es.clusters <- Idents(se.integrated)
  
  dir.create('escape')
  geneset.c5 <- getGeneSets(species = species, library = 'C5', subcategory = 'BP')

  # as a heuristic and for speed, removing genesets >1500 
  size <- NULL
  for (x in seq_along(geneset.c5)) {
    size <- c(size, length(geneset.c5[[x]]))
  }
  remove <- unname(which(size >= 1500))
  if (length(remove) != 0) {
    geneset.c5 <- geneset.c5[-remove]
  }
  
  if (length(pathways) == 1 & pathways[1] == "none"){
    for (i in 1:length(pathways)){
      if (any(grepl(pathways[i], names(geneset.c5)))){
        print(paste0(pathways[i], " exists in gene set, proceed"))
      } else {
        stop(paste0("no matching pathways for ", pathways[i]))
      }
    }
  }
  
  print("Start scGSEA analysis, this will take a while")
  se.integrated <- runEscape(se.integrated, 
                             method = "UCell",
                             gene.sets = geneset.c5, 
                             normalize = FALSE,
                             groups = 5000,
                             new.assay.name = "escape.UCell")

  print("Normalizing UCell values")
  se.integrated <- performNormalization(se.integrated, 
                                        assay = "escape.UCell", 
                                        gene.sets = geneset.c5, 
                                        scale.factor = se.integrated$nFeature_RNA)
  saveRDS(se.integrated, "se_integrated_escape_norm.rds")

  print("Performing DE analysis for the pathways")
  ucell.markers <- FindAllMarkers(se.integrated, assay = "escape.UCell_normalized")

  if (length(pathways) == 1 & pathways[1] == "none"){
    ucell.markers.top5 <-  as.data.frame(ucell.markers %>% filter(p_val_adj <= 0.05) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
    if (dim(ucell.markers.top5)[1] != 0 & dim(ucell.markers.top5)[2] != 0){
      if (length(unique(ucell.markers.top5$gene)) > 1){
        p1 <- heatmapEnrichment(se.integrated,
                                group.by = "es.clusters",
                                gene.set.use = unique(ucell.markers.top5$gene),
                                assay = "escape.UCell_normalized",
                                scale = TRUE) +
          theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
          ggtitle("Top 5 GO Pathways per Cluster")

        ggSaveAndJPEG(p1, 
                     paste0("escape/", "escape_heatmap_top5"),
                     plots.format = plots.format,
                     width = 20,
                     height = 10
        )
        
      } else {
        print("Only 1 significant pathway")
      }
      
      for (i in 1:nrow(ucell.markers.top5)){
        p2 <- geyserEnrichment(se.integrated, 
                               assay = "escape.UCell_normalized",
                               group.by = "es.clusters",
                               gene.set = ucell.markers.top5[i,7], 
                               facet.by = "group",
                               scale = TRUE) +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
          ggtitle(paste0("Geyser Enrichment: ", levels(droplevels(ucell.markers.top5[i,6])))) + 
          theme(legend.position = "none")
        tar.dir <- paste0("escape/", levels(droplevels(ucell.markers.top5[i,6])), "/")
        dir.create(tar.dir)
        
        ggSaveAndJPEG(p2, 
                     paste0(tar.dir, unique(ucell.markers.specific.filtered$gene)[j], "_geyser"),
                     plots.format = plots.format,
                     width = 12,
                     height = 12
        )
      }
    } else {
      print("No significant pathways")
    }
  } else {
    ## if specified pathways were wanted, folders will be split by the general phrase and then filled with geyser plots
    for (i in 1:length(pathways)){
      if (any(grepl(pathways[i], ucell.markers$gene))){
        ucell.markers.specific <- ucell.markers[grep(pathways[i], ucell.markers$gene),]
        ucell.markers.specific.filtered <- as.data.frame(ucell.markers.specific %>% filter(p_val_adj <= 0.05) %>% group_by(cluster))
        
        if (dim(ucell.markers.specific.filtered)[1] != 0 & dim(ucell.markers.specific.filtered)[2] != 0){
          if (length(unique(ucell.markers.specific.filtered$gene)) > 1){
            p1 <- heatmapEnrichment(se.integrated,
                                    group.by = "es.clusters",
                                    gene.set.use = unique(ucell.markers.specific.filtered$gene),
                                    assay = "escape.UCell_normalized",
                                    scale = TRUE) +
              theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
              ggtitle(paste0(pathways[i], " GO related Pathways"))
            
            ggSaveAndJPEG(p1, 
                         paste0("escape/", "escape_heatmap_", pathways[i]),
                         plots.format = plots.format,
                         width = 20,
                         height = 10
                         )
            
          } else {
            print("Only 1 significant pathway")
          }
          
          for (j in 1:length(unique(ucell.markers.specific.filtered$gene))){
            p2 <- geyserEnrichment(se.integrated, 
                                   assay = "escape.UCell_normalized",
                                   group.by = "es.clusters",
                                   gene.set = unique(ucell.markers.specific.filtered$gene)[j], 
                                   facet.by = "group",
                                   scale = TRUE) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
              ggtitle(paste0("Geyser Enrichment")) + 
              theme(legend.position = "none")
            tar.dir <- paste0("escape/", pathways[i], "/")
            dir.create(tar.dir)
            ggSaveAndJPEG(p2, 
                         paste0(tar.dir, unique(ucell.markers.specific.filtered$gene)[j], "_geyser"),
                         plots.format = plots.format,
                         width = 12,
                         height = 12
                         )
          }
        } else {
          print(paste0("No significant pathways for "), pathways[i])
        }
      } else {
        print("PATHWAY WITH PHRASE DOES NOT EXIST: rerun the escape analyses with with proper name")
      }
    }
  }
  se.integrated$es.clusters <- NULL
}
