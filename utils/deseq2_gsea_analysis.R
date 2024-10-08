set.seed(333)

DESeq2ConditionPerCluster <-  function(se.integrated, species){
  se.integrated$de.clusters <- Idents(se.integrated)
  if (length(se.integrated@misc) != 0 ){
    bulk <- AggregateExpression(se.integrated, return.seurat = T, 
                                assays = "RNA", 
                                group.by = c("de.clusters", "sample", "group", se.integrated@misc[["co.conditions"]]))
    #group.by = c("de.clusters", "sample", "group", se.integrated@misc[["co.conditions"]]))
  } else {
    bulk <- AggregateExpression(se.integrated, return.seurat = T, 
                                assays = "RNA", 
                                group.by = c("de.clusters", "sample", "group"))
  }
  #bulk$de.clusters <- Idents(bulk)
  Idents(bulk) <-   bulk$de.clusters
  
  m_df<- msigdbr(species = species, category = "C5", subcategory = "BP") # dont need to reload the dataset every time for GSEA
  fgsea_sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  for (i in 1:nlevels(se.integrated@meta.data[["de.clusters"]])){
    #print(paste, 'Performing DE analysis on cluster:', i)
    cluster.bulk <- NA
    cluster.name <- NA
    if (all(as.character(se.integrated$de.clusters) == as.character(se.integrated$seurat_clusters))){
      cluster.name <- paste0("g", i-1)
      print(cluster.name)
      cluster.bulk <- subset(bulk, de.clusters == cluster.name) # seurat subset doesnt seem to like string args

    } else {
      cluster.name <- levels(droplevels(se.integrated@meta.data[["de.clusters"]]))[i]
      print(cluster.name)
    	if(grepl("_", cluster.name, fixed=TRUE)){
    		cluster.name <- sub("_", "-", cluster.name)
    	}
      cluster.bulk <- subset(bulk, de.clusters == cluster.name)

    }
    
    group.pairs <- NA
    list.comparisons <- list("group")
    if (length(se.integrated@misc) != 0 ){
      for (j in se.integrated@misc$co.conditions){
        group.co <- paste0("group_", j)
        cluster.bulk <- AddMetaData(cluster.bulk, paste(cluster.bulk$group, cluster.bulk@meta.data[[j]], sep = "_"), group.co)
        list.comparisons <- append(list.comparisons, group.co)
      }
    }
    
    for (k in 1:length(list.comparisons)){
      Idents(cluster.bulk) <- list.comparisons[[k]]
      grouping <- list.comparisons[[k]]
      if (any(duplicated(as.data.frame(cluster.bulk@meta.data[[grouping]])))){ # check if you have at least more than 1 sample for the comparison
        
        # create pairs for each possible combination
        group.pairs <- as.data.frame(combn(unique(cluster.bulk@meta.data[[grouping]]), 2))
        group.pairs <- sapply(group.pairs, function(x) as.character(x), simplify = FALSE)
        
        if (k >= 2){
          group.pairs.filt <- list()
          for (z in group.pairs){
            time.match <- sub("^[^_]*_", "", z)
            if (time.match[1] == time.match[2]) {
              cond.mismatch <- sub("_.*", "", z)
              if (cond.mismatch[1] != cond.mismatch[2]){
                group.pairs.filt <- append(group.pairs.filt, list(z))
              }
            }
          }
          group.pairs <- group.pairs.filt
          print(group.pairs)
        }
        
        for (j in 1:length(group.pairs)){ # perform each pairwise comparison
          target <- group.pairs[[j]]
          print(target[1])
          print(target[2])
          
          dup_check <- as.data.frame(cluster.bulk@active.ident)
          colnames(dup_check)[1] <- "cond"
          value_counts <- table(as.character(dup_check$cond))
          
          # Extract values that occur only once
          unique_values <- names(value_counts[value_counts == 1])
          
          if (length(unique_values) == 0){
            de_markers <- FindMarkers(cluster.bulk, ident.1 = target[1], ident.2 = target[2], slot = "counts", test.use = "DESeq2",
                                      verbose = T, min.cells.feature = 0, min.cells.group = 0)
            # start GSEA analysis here too since will be doing all the comparisons here
            de_markers$gene <- rownames(de_markers)
            GseaComparison(de_markers, cluster.name, target[1], target[2], fgsea_sets)
            
            #print('plot')
            p <- ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + 
              geom_point(size = 0.5, alpha = 0.5) + 
              theme_bw() +
              ylab("-log10(unadjusted p-value)") + 
              geom_text_repel(aes(label = ifelse(((p_val_adj < 0.05 & avg_log2FC >= 2)|(p_val_adj < 0.05 & avg_log2FC <= -2)), gene,
                                                                                      "")), colour = "red", size = 3) + 
              ggtitle(paste0("DESeq2: ", cluster.name, " ", target[1], " vs ", target[2]))
            #dir.create(paste('deseq2', sep=''))
            #deseq2.folder <- paste(plot.path, 'deseq2/')
            #print('save')
            write.table(de_markers, paste0("deseq2_cluster_", cluster.name, "_", target[1], "_vs_", target[2], '.txt'), row.names = F, quote = F)
            PrintSave(p, paste0("deseq2_cluster_", cluster.name, "_", target[1], "_vs_", target[2], '.pdf'))
          }
        }
      } else {
        print(paste0("No comparison analysis available for ", k, " since you only have 1 replicate, skipping this"))
      }
    }
  }
  se.integrated$de.clusters <- NULL
}

GseaComparison <- function(de.markers, cluster.name, ident.1, ident.2, fgsea.sets){
  cluster.genes <- de.markers %>%
    arrange(desc(avg_log2FC)) %>% 
    dplyr::select(gene, avg_log2FC) # use avg_log2FC as ranking for now; https://www.biostars.org/p/9526168/
  
  ranks <- deframe(cluster.genes)
  
  fgseaRes <- fgsea(fgsea.sets, stats = ranks)
  fgseaRes$leadingEdge <- as.character(fgseaRes$leadingEdge)
  write.table(fgseaRes, paste0("gsea_cluster_", cluster.name, "_", ident.1, "_vs_", ident.2,".txt"), quote = FALSE,row.names = T, sep = "\t", col.names = T)
  fgseaRes <- filter(fgseaRes, padj <= 0.05 & size >= 3) %>% arrange(desc(NES))
  fgseaRes$Enrichment = ifelse(fgseaRes$NES > 0, "red", "blue") 
  
  filtRes <-  rbind(head(fgseaRes, n = 10),
                    tail(fgseaRes, n = 10 ))
  
  p1 <- ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point(aes(fill = Enrichment, size = size), shape=21) + # size is equal to the amount of genes found in the given gene set
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    theme_bw() +
    labs(x="Pathway", y="Normalized Enrichment Score") + 
    scale_fill_identity() + 
    theme(legend.position = 'none') + 
    ggtitle(paste0("GSEA: ", cluster.name, " ", ident.1, " vs ", ident.2))
  
  #dir.create(paste(plot.path, 'gsea', sep=''))
  PrintSave(p1, paste0("gsea_cluster_", cluster.name, "_", ident.1, "_vs_", ident.2, '.pdf'), w=12)
}
