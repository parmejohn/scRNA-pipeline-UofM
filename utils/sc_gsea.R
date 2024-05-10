set.seed(333)

EscapeGSEA <- function(se.integrated, species, pathways){
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
      if (grepl(pathways[i], names(geneset.c5))){
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
    
    p1 <- heatmapEnrichment(se.integrated,
                            group.by = "ident",
                            gene.set.use = unique(ucell.markers.top5$gene),
                            assay = "escape.UCell_normalized",
                            scale = TRUE) +
      theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
      ggtitle("Top 5 GO Pathways per Cluster")
    ggsave(paste0("escape/", "escape_heatmap_top5.pdf"), p1, width = 20, height = 10)
    
    for (i in 1:nrow(ucell.markers.top5)){
      p2 <- geyserEnrichment(se.integrated, 
                             assay = "escape.UCell_normalized",
                             group.by = "ident",
                             gene.set = ucell.markers.top5[i,7], 
                             facet.by = "group",
                             scale = TRUE) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        ggtitle(paste0("Geyser Enrichment: ", levels(droplevels(ucell.markers.top5[i,6])))) + 
        theme(legend.position = "none")
      tar.dir <- paste0("escape/", levels(droplevels(ucell.markers.top5[i,6])), "/")
      dir.create(tar.dir)
      ggsave(paste0(tar.dir, ucell.markers.top5[i,7], "_geyser.pdf"), p2, width = 12, height = 12)
    }
  } else {
    ## if specified pathways were wanted, folders will be split by the general phrase and then filled with geyser plots
    for (i in 1:length(pathways)){
      if (grepl(pathways[i], ucell.markers$gene)){
        ucell.markers.specific <- ucell.markers[grep(pathways[i], ucell.markers$gene),]
        ucell.markers.specific.filtered <- as.data.frame(ucell.markers.specific %>% filter(p_val_adj <= 0.05) %>% group_by(cluster))
        
        p1 <- heatmapEnrichment(se.integrated,
                                group.by = "ident",
                                gene.set.use = unique(ucell.markers.specific.filtered$gene),
                                assay = "escape.UCell_normalized",
                                scale = TRUE) +
          theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
          ggtitle(paste0(pathways[i], " GO related Pathways"))
        ggsave(paste0("escape/", "escape_heatmap_", pathways[i], ".pdf"), p1, width = 20, height = 10)
        
        for (j in 1:length(unique(ucell.markers.specific.filtered$gene))){
          p2 <- geyserEnrichment(se.integrated, 
                                 assay = "escape.UCell_normalized",
                                 group.by = "ident",
                                 gene.set = unique(ucell.markers.specific.filtered$gene)[j], 
                                 facet.by = "group",
                                 scale = TRUE) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
            ggtitle(paste0("Geyser Enrichment")) + 
            theme(legend.position = "none")
          tar.dir <- paste0("escape/", pathways[i], "/")
          dir.create(tar.dir)
          ggsave(paste0(tar.dir, unique(ucell.markers.specific.filtered$gene)[j], 
                        "_geyser.pdf"), p2, width = 12, height = 12)
        }
      } else {
        print("PATHWAY WITH PHRASE DOES NOT EXIST: rerun the escape analyses with with proper name")
      }
    }
  }
}
