#source(paste0(dirname(dirname(dirname(getwd()))),"/utils/misc.R"))
set.seed(333)


DESeq2ConditionPerCluster <-  function(se.integrated, species){
  se.integrated$de.clusters <- Idents(se.integrated)
  bulk <- AggregateExpression(se.integrated, return.seurat = T, 
                              assays = "RNA", 
                              group.by = c("de.clusters", "sample", "group"))
  #bulk$de.clusters <- Idents(bulk)
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
      #expr <- FetchData(bulk, vars = "de.clusters")
      #cluster.bulk <- bulk[, which(expr == cluster.name)]
    } else {
      cluster.name <- levels(droplevels(se.integrated@meta.data[["de.clusters"]]))[i] ### NEED TO SEE IF THIS WORKS
      print(cluster.name)
	if(grepl("_", cluster.name, fixed=TRUE)){
		cluster.name <- sub("_", "-", cluster.name)
		print(cluster.name)
	}
      cluster.bulk <- subset(bulk, de.clusters == cluster.name)
      #expr <- FetchData(bulk, vars = "de.clusters")
      #cluster.bulk <- bulk[, which(expr == cluster.name)]
    }
    
    Idents(cluster.bulk) <- "group"
    
    # create pairs for each possible combination
    group.pairs <- as.data.frame(combn(unique(se.integrated@meta.data[["group"]]), 2))
    group.pairs <- sapply(group.pairs, function(x) as.character(x), simplify = FALSE)
    
    for (j in 1:length(group.pairs)){ # perform each pairwise comparison
      target <- group.pairs[[j]]
      de_markers <- FindMarkers(cluster.bulk, ident.1 = target[1], ident.2 = target[2], slot = "counts", test.use = "DESeq2",
                                verbose = T, min.cells.feature = 0, min.cells.group = 0)
      # start GSEA analysis here too since will be doing all the comparisons here
      de_markers$gene <- rownames(de_markers)
      GseaComparison(de_markers, cluster.name, target[1], target[2], fgsea_sets)
      
      #print('plot')
      p <- ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
        ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(((p_val_adj < 0.05 & avg_log2FC >= 2)|(p_val_adj < 0.05 & avg_log2FC <= -2)), gene,
                                                                                "")), colour = "red", size = 3)
      #dir.create(paste('deseq2', sep=''))
      #deseq2.folder <- paste(plot.path, 'deseq2/')
      #print('save')
      write.table(de_markers, paste0("deseq2_cluster_", cluster.name, "_", target[1], "_vs_", target[2], '.txt'), row.names = F, quote = F)
      PrintSave(p, paste0("deseq2_cluster_", cluster.name, "_", target[1], "_vs_", target[2], '.pdf'))
    }
  }
  se.integrated$de.clusters <- NULL
  
}

GseaComparison <- function(de.markers, cluster.name, ident.1, ident.2, fgsea.sets){
  cluster.genes<- de.markers %>%
    arrange(desc(avg_log2FC)) %>% 
    dplyr::select(gene, avg_log2FC) # use avg_log2FC as ranking for now; https://www.biostars.org/p/9526168/
  
  ranks <- deframe(cluster.genes)
  
  fgseaRes <- fgsea(fgsea.sets, stats = ranks, nperm = 1000)
  fgseaRes <- filter(fgseaRes, pval <= 0.01) %>% arrange(desc(NES))
  fgseaRes$Enrichment = ifelse(fgseaRes$NES > 0, "Up-regulated", "Down-regulated") 
  
  filtRes <-  rbind(head(fgseaRes, n = 10),
                    tail(fgseaRes, n = 10 ))
  
  p1 <- ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) + # size is equal to the amount of genes found in the given gene set
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score")
  
  #dir.create(paste(plot.path, 'gsea', sep=''))
  PrintSave(p1, paste0("gsea_cluster_", cluster.name, "_", ident.1, "_vs_", ident.2, '.pdf'), w=12)
}
