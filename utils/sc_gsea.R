#source(paste0(dirname(dirname(dirname(getwd()))),"/utils/misc.R"))
set.seed(333)

EscapeGSEA <- function(se.integrated, species){
  dir.create('escape')
  geneset.c5 <- getGeneSets(species = species, library = 'C5', subcategory = 'BP')
  #geneset.c5.bp <- getGeneSets(species = 'Mus musculus', library = 'C5', subcategory = "BP")
  
  # as a heuristic and for speed, removing genesets >1500 
  size <- NULL
  for (x in seq_along(geneset.c5)) {
    size <- c(size, length(geneset.c5[[x]]))
  }
  remove <- unname(which(size >= 1500))
  if (length(remove) != 0) {
    geneset.c5 <- geneset.c5[-remove]
  }
  
  print("Start scGSEA analysis, this will take a while")
  se.integrated <- runEscape(se.integrated, 
                             method = "UCell",
                             gene.sets = geneset.c5, 
                             normalize = FALSE,
                             groups = 5000,
                             new.assay.name = "escape.UCell")
  saveRDS(se.integrated, "se_integrated_escape.rds")
  
  # takes over 8+ hours
  print("Normalizing UCell values")
  se.integrated <- performNormalization(se.integrated, 
                                        assay = "escape.UCell", 
                                        gene.sets = geneset.c5, 
                                        scale.factor = se.integrated$nFeature_RNA)
  saveRDS(se.integrated, "se_integrated_escape_norm.rds")
  #print(paste0("Saved object under ", out.path, "/se_integrated_escape_norm.rds"))
  
  print("Performing DE analysis for the pathways")
  ucell.markers <- FindAllMarkers(se.integrated, assay = "escape.UCell_normalized")
  ucell.markers.top5 <-  as.data.frame(ucell.markers %>% filter(p_val_adj <= 0.05) %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC))
  
  p1 <- heatmapEnrichment(se.integrated,
                          group.by = "ident",
                          gene.set.use = unique(ucell.markers.top5$gene),
                          assay = "escape.UCell_normalized",
                          scale = TRUE) +
    theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(paste0("escape/", "escape_heatmap_top5.pdf"), p1, width = 20, height = 10)
  
  for (i in 1:nrow(ucell.markers.top5)){
    p2 <- geyserEnrichment(se.integrated, 
                           assay = "escape.UCell_normalized",
                           group.by = "ident",
                           gene.set = ucell.markers.top5[i,7], 
                           facet.by = "group",
                           scale = TRUE) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    tar.dir <- paste0("escape/", levels(droplevels(ucell.markers.top5[i,6])), "/")
    dir.create(tar.dir)
    ggsave(paste0(tar.dir, ucell.markers.top5[i,7], "_geyser.pdf"), p2, width = 12, height = 12)
  }
}
