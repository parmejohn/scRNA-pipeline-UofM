AtacAnalyses <- function(se.integrated.atac, species){
  m_df<- msigdbr(species = species, category = "C5", subcategory = "BP") # dont need to reload the dataset every time for GSEA
  fgsea.sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  DefaultAssay(se.integrated.atac) <- "ATAC"
  
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
  )
  
  if (species == "Homo sapiens"){
    genome = BSgenome.Hsapiens.UCSC.hg38
  } else if (species == "Mus musculus"){
    genome = BSgenome.Mmusculus.UCSC.mm10
  } else {
    stop("invalid species")
  }
  
  se.integrated.atac <- AddMotifs(
    object = se.integrated.atac,
    genome = genome,
    pfm = pfm
  )
  
  se.integrated.atac$dap_clusters <- Idents(se.integrated.atac)
  Idents(se.integrated.atac) <- "group"
  
  se.integrated.atac <- RegionStats(se.integrated.atac, genome = BSgenome.Hsapiens.UCSC.hg38)
  
  dir.create("dap_plots")
  dir.create("dap_data")
  
  for (i in 1:nlevels(se.integrated.atac@meta.data[["dap_clusters"]])){
    
    if (all(as.character(se.integrated.atac$dap_clusters) == as.character(se.integrated.atac$seurat_clusters))){
      cluster.name <- i - 1
      se.integrated.atac.filt <- subset(se.integrated.atac, dap_clusters == cluster.name) # seurat subset doesnt seem to like string args
      
    } else {
      cluster.name <- levels(droplevels(se.integrated.atac@meta.data[["dap_clusters"]]))[i]
      se.integrated.atac.filt <- subset(se.integrated.atac, dap_clusters == cluster.name)
    }
    
    print(paste0("on cluster ", cluster.name))
    
    list.comparisons <- ListAllPossibleComparisons(se.integrated.atac = se.integrated.atac,
                                                   se.integrated.atac.filt = se.integrated.atac.filt)
    for (k in 1:length(list.comparisons)){
      Idents(se.integrated.atac.filt) <- list.comparisons[[k]]
      grouping <- list.comparisons[[k]]
      if (any(duplicated(as.data.frame(se.integrated.atac.filt@meta.data[[grouping]])))){ # check if you have at least more than 1 sample for the comparison
        group.pairs <- MatchCovariantGroupPairs(se.integrated.atac.filt = se.integrated.atac.filt,
                                                grouping = grouping,
                                                not.main.group = k)
        
        for (j in 1:length(group.pairs)){ # perform each pairwise comparison
          target <- group.pairs[[j]]
          print(target[1])
          print(target[2])
          
          data.dir <- paste0("dap_data/", target[1], "_vs_", target[2], "/")
          dir.create(data.dir)
          
          plots.dir <- paste0("dap_plots/", target[1], "_vs_", target[2], "/")
          dir.create(plots.dir)
          
          volcano_plots <- paste0(plots.dir, "volcano_plots/")
          dir.create(volcano_plots)
          
          volcano_data <- paste0(data.dir, "volcano_data/")
          dir.create(volcano_data)
          
          closest_gene_plots <- paste0(plots.dir, "closest_gene_plots/")
          dir.create(closest_gene_plots)
          
          motif_plots <- paste0(plots.dir, "motif_plots/")
          dir.create(motif_plots)
          
          da.peaks.all <- dap_volcano_plot(se.integrated.atac.filt, target[1], target[2], cluster.name = cluster.name, plots.dir = volcano_plots, data.dir = volcano_data)
          print("volcano")
          
          ### DAP specific analyses
          # pulling top and bot 5 of DAPs, also making sure that the are up or dn regulated
          # plot out gene expression
          closest.gene.dap.up <- FindTopDAPGenes(da.peaks.all = da.peaks.all, se.integrated.atac = se.integrated.atac.filt,
                                                 cluster.name = cluster.name, group = target[1], top = 5, plots.dir = closest_gene_plots)
          closest.gene.dap.dn <- FindTopDAPGenes(da.peaks.all = da.peaks.all, se.integrated.atac = se.integrated.atac.filt,
                                                 cluster.name = cluster.name, group = target[2], top = -5, plots.dir = closest_gene_plots)
          print("expression")
          
          
          # perform GSEA, ranked by log2fc 
          all.closest <- ClosestFeature(se.integrated.atac.filt, regions = da.peaks.all$gene)
          colnames(all.closest)[7] <- "gene"
          all.closest <- merge(all.closest, da.peaks.all, by = "gene")
          
          GseaAtac(all.closest, fgsea.sets, cluster.name = cluster.name, ident.1 = target[1], ident.2 = target[2], plots.dir = closest_gene_plots) 
          print("gsea")
          
          ## coverage plot of the DAPs
          CoveragePlotDAP(se.integrated.atac = se.integrated.atac.filt, closest.gene.dap = closest.gene.dap.up, cluster.name = cluster.name, plots.dir = closest_gene_plots)
          CoveragePlotDAP(se.integrated.atac = se.integrated.atac.filt, closest.gene.dap = closest.gene.dap.dn, cluster.name = cluster.name, plots.dir = closest_gene_plots)
          print("covearage")
          
          
          ### DNA motifs that are overrepresented in a set of peaks that are differentially accessible between cell types
          enriched.motifs.up <- find_enriched_motifs(se.integrated.atac.filt = se.integrated.atac.filt, 
                                                     da.peaks.all = da.peaks.all, cutoff = 0.2, pos = T, group = target[1], 
                                                     cluster.name = cluster.name, plots.dir = motif_plots)
          enriched.motifs.down <- find_enriched_motifs(se.integrated.atac.filt = se.integrated.atac.filt, 
                                                       da.peaks.all = da.peaks.all, cutoff = 0.2, pos = F, group = target[2], 
                                                       cluster.name = cluster.name, plots.dir = motif_plots)
          print("motifs")
          
          
          # gather the footprinting information for sets of motifs; leaving out for now because of computational constraints
          # TopMotifFootprints(se.integrated.atac = se.integrated.atac.filt, enriched.motifs = enriched.motifs.up, cluster.name = i, group = target[1], plots.dir = motif_plots)
          # TopMotifFootprints(se.integrated.atac = se.integrated.atac.filt, enriched.motifs = enriched.motifs.down, cluster.name = i, group = target[2], plots.dir = motif_plots)
          print("fp")
        }
      } else {
        print(paste0("No comparison analysis available for ", k, " since you only have 1 replicate, skipping this"))
      } 
    }
  }
}

dap_volcano_plot <- function(se.object, ident.1, ident.2, test.use = "wilcox", min.pct = 0.05, cluster.name, p_val_adj_cutoff = 0.05, avg_log2FC_cutoff = 1.5, plots.dir, data.dir){
  da.peaks <- RunPresto(
    object = se.object,
    ident.1 = ident.1,
    ident.2 = ident.2,
    min.pct = min.pct,
    logfc.threshold = 0.1
  )
  da.peaks$gene <- rownames(da.peaks)
  write.table(da.peaks, paste0(data.dir ,"dap_", "cluster_", cluster.name, ident.1, "_vs_", ident.2, ".txt"), 
              quote = FALSE, row.names = T, sep = "\t", col.names = T)
  
  da.peaks$sig <- "Not Significant"
  da.peaks$sig[da.peaks[["avg_log2FC"]] >= avg_log2FC_cutoff & 
                 da.peaks[["p_val_adj"]] < p_val_adj_cutoff] <- "Upregulated"
  da.peaks$sig[da.peaks[["avg_log2FC"]] < -(avg_log2FC_cutoff) & 
                 da.peaks[["p_val_adj"]] <= p_val_adj_cutoff] <- "Downregulated"
  da.peaks$sig <- factor(da.peaks$sig, levels = c("Downregulated", "Not Significant", "Upregulated"))
  
  p <- ggplot(da.peaks, aes(avg_log2FC, -log10(p_val), color = sig)) + 
    geom_point(size = 0.5, alpha = 0.5) + 
    scale_color_manual(values = c("blue", "gray", "red")) +
    theme_bw() +
    ylab("-log10(unadjusted p-value)") + 
    ggtitle(paste0("DAPs: ", "cluster ", cluster.name, " ", ident.1, " vs ", ident.2))
  ggsave(paste0(plots.dir, "scatac_volcano_", "cluster_", cluster.name, "_", ident.1, "_vs_", ident.2, ".pdf"), plot = p, width = 8, height = 8)
  
  return(da.peaks)
}

GseaAtac <- function(all.closest, fgsea.sets, padj.cutoff = 0.05, cluster.name, ident.1, ident.2, plots.dir) {
  cluster.genes <- all.closest %>%
    arrange(desc(avg_log2FC)) %>% 
    dplyr::select(gene_name, avg_log2FC) # use avg_log2FC as ranking for now; https://www.biostars.org/p/9526168/
  
  ranks <- deframe(cluster.genes)
  
  fgseaRes <- fgsea(fgsea.sets, stats = ranks)
  fgseaRes$leadingEdge <- as.character(fgseaRes$leadingEdge)
  
  fgseaRes <- fgsea(fgsea.sets, stats = ranks)
  fgseaRes$leadingEdge <- as.character(fgseaRes$leadingEdge)
  
  fgseaRes <- filter(fgseaRes, padj < padj.cutoff & size >= 3) %>% arrange(desc(NES))
  fgseaRes$Enrichment = ifelse(fgseaRes$NES > 0, "blue", "red") 
  
  filtRes <-  rbind(head(fgseaRes, n = 10),
                    tail(fgseaRes, n = 10 ))
  
  p <- ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point(aes(fill = Enrichment, size = size), shape=21) + # size is equal to the amount of genes found in the given gene set
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    theme_bw() +
    labs(x="Pathway", y="Normalized Enrichment Score") + 
    scale_fill_identity() + 
    theme(legend.position = 'none') + 
    ggtitle(paste0("GSEA: ", cluster.name, " ", ident.1, " vs ", ident.2))
  ggsave(paste0(plots.dir, "scatac_closest_genes_dap_gsea_", "cluster_", cluster.name, ".pdf"), plot = p, width = 8, height = 8)
}

FindTopDAPGenes <- function(da.peaks.all, p_val_adj.cutoff = 0.05, avg_log2FC.cutoff = 1, se.integrated.atac, cluster.name, group, top, plots.dir) {
  da.peaks.top  <- as.data.frame(da.peaks.all %>%
                                   subset(p_val_adj < p_val_adj.cutoff & abs(avg_log2FC) >= avg_log2FC.cutoff) %>%
                                   top_n(top, avg_log2FC))
  if (top > 0){
    da.peaks.top <- filter(da.peaks.top, avg_log2FC > 0)
  } else {
    da.peaks.top <- filter(da.peaks.top, avg_log2FC < 0)
  }
  if (nrow(da.peaks.top) > 0) {
    closest.gene.dap <- ClosestFeature(se.integrated.atac, regions = da.peaks.top$gene)
    print(closest.gene.dap)
    print("filter genes")
    for (v in 1:nrow(closest.gene.dap)){
      if(!(closest.gene.dap$gene_name[v] %in% rownames(se.integrated.atac@assays[["RNA"]]@features) & 
           closest.gene.dap$gene_name[v] %in% Annotation(se.integrated.atac)$gene_name)){
        closest.gene.dap <- closest.gene.dap[-v,]
      }
    }
    
    
    if (nrow(closest.gene.dap) > 0){
      p <- VlnPlot(
        object = se.integrated.atac,
        features = closest.gene.dap$gene_name,
        pt.size = 0.1
      ) + plot_annotation(paste0("Cluster ", cluster.name, " ", group))
      ggsave(paste0(plots.dir, "scatac_closest_genes_dap_gex_", "cluster_", cluster.name, "_", group, ".pdf"), plot = p, width = 8, height = 8)
    }
  } else {
    closest.gene.dap <- data.frame()
  }
  return(closest.gene.dap)
}

CoveragePlotDAP <- function(se.integrated.atac, closest.gene.dap, cluster.name, plots.dir, min.cells = 10) {
  # calculate the links between gene and peaks
  print(closest.gene.dap)
  if (nrow(closest.gene.dap) > 0){
    # if (any(grepl("RP*.*", closest.gene.dap$gene_name))){
    #   closest.gene.dap <- closest.gene.dap[- grep("RP*.*", closest.gene.dap$gene_name),]
    # }
    
    print("linking")
    
    expression.data <- GetAssayData(
      object = se.integrated.atac, assay = "RNA", layer = "data"
    )
    genecounts <- rowSums(x = expression.data > 0)
    genes.keep <- genecounts > min.cells
    genes.keep <- intersect(
      x = names(x = genes.keep[genes.keep]), y = closest.gene.dap$gene_name
    )

	# checking if peaks are in distance like in the LinkPeaks function
	peak.data <- GetAssayData(
    		object = se.integrated.atac, assay = "ATAC", layer = "counts"
  	)
	peakcounts <- rowSums(x = peak.data > 0)
	peaks.keep <- peakcounts > min.cells
	peak.data <- peak.data[peaks.keep, ]

	peaks <- granges(x = se.integrated.atac[["ATAC"]])
  	peaks <- peaks[peaks.keep]
    	annot <- Annotation(object = se.integrated.atac[["ATAC"]])
   	if (is.null(x = annot)) {
      		stop("Gene annotations not found")
    	}
    	gene.coords <- CollapseToLongestTranscriptInternal(
      		ranges = annot
    	)
	expression.data <- expression.data[genes.keep, , drop = FALSE]
	genes <- rownames(x = expression.data)
	gene.coords.use <- gene.coords[gene.coords$gene_name %in% genes,]


    peak_distance_matrix <- DistanceToTSSInternal(
    	peaks = peaks,
    	genes = gene.coords.use,
    	distance = 5e+05
   )

    if (length(genes.keep) > 0 && sum(peak_distance_matrix) != 0){
      se.integrated.atac <- LinkPeaks(
        object = se.integrated.atac,
        peak.assay = "ATAC",
        expression.assay = "RNA",
        genes.use = c(genes.keep, "PRDM1"),
        min.cells = min.cells
      )
    }
    
    # plot out a coverage plot for the top 5 and bottom 5 genes
    # highlight the DAP region if it is found within the gene region
    for (v in 1:nrow(closest.gene.dap)){
      #          if (closest.gene.dap.up$distance[v] == 0){ # if the highlight is not in the region it wont show; no need for if statement
      print(closest.gene.dap$gene_name[v])
      regions_highlight <- subsetByOverlaps(StringToGRanges(closest.gene.dap$query_region[v]), 
                                            LookupGeneCoords(se.integrated.atac, closest.gene.dap$gene_name[v]))
      print("plot cover")
      if(closest.gene.dap$gene_name[v] %in% rownames(se.integrated.atac@assays[["RNA"]]@features) & 
         closest.gene.dap$gene_name[v] %in% Annotation(se.integrated.atac)$gene_name){
        p <- CoveragePlot(
          object = se.integrated.atac,
          region = closest.gene.dap$gene_name[v],
          features = closest.gene.dap$gene_name[v],
          extend.upstream = 1000,
          extend.downstream = 1000,
          assay = 'ATAC', 
          expression.assay = 'RNA',
          region.highlight = regions_highlight
        ) + plot_annotation(paste0("Cluster ", cluster.name, ", DAP: ", closest.gene.dap$query_region[v]))
        print("plot done")
        ggsave(paste0(plots.dir, "scatac_closest_genes_dap_coverage_", "cluster_", cluster.name, "_", closest.gene.dap$gene_name[v],".pdf"), 
               plot = p,
               width = 8, 
               height = 8)
        #         }
      }
    }
  }
}

find_enriched_motifs <- function(se.integrated.atac.filt, da.peaks.all, p_val_adj_cutoff = 0.05, cutoff = 0.2, pos, group, cluster.name, plots.dir){
  if (pos){
    da.peaks.sub <- subset(da.peaks.all, p_val_adj < p_val_adj_cutoff & 
                             pct.1 > cutoff & 
                             avg_log2FC > 0
    )
  } else {
    da.peaks.sub <- subset(da.peaks.all, p_val_adj < p_val_adj_cutoff & 
                             pct.2 > cutoff &
                             avg_log2FC < 0
    )
  }
  if (any(grepl("^KI", da.peaks.sub$gene))){
    da.peaks.sub <- da.peaks.sub[- grep("KI", da.peaks.sub$gene),]
  }
  if (any(grepl("^GL", da.peaks.sub$gene))){
    da.peaks.sub <- da.peaks.sub[- grep("GL", da.peaks.sub$gene),]
  }
  print(da.peaks.sub)
  if(nrow(da.peaks.sub) > 0){
    enriched.motifs.sub <- FindMotifs(
      object = se.integrated.atac.filt,
      features = da.peaks.sub$gene
    )
    p <- MotifPlot(
      object = se.integrated.atac.filt,
      motifs = head(rownames(enriched.motifs.sub), 6)
    )
    ggsave(paste0(plots.dir, "scatac_motif_", "cluster_", cluster.name, "_", group, ".pdf"), plot = p, width = 8, height = 8)
  } else {
    enriched.motifs.sub <- data.frame()
  }
  return(enriched.motifs.sub)
}

TopMotifFootprints <- function(se.integrated.atac, enriched.motifs, cluster.name, group, plots.dir) {
  if (nrow(enriched.motifs) > 0){
    se.integrated.atac <- Footprint(
      object = se.integrated.atac,
      motif.name = head(rownames(enriched.motifs), 6),
      genome = BSgenome.Hsapiens.UCSC.hg38,
      in.peaks = TRUE,
      compute.expected = FALSE
    )
    
    #saveRDS(se.integrated.atac.filt.fp, file = "/research/2024_scatac/rscript_ex_outputs/se_integrated_atac_fp.rds")
    
    # plot the footprint data for each group of cells
    p <- PlotFootprint(se.integrated.atac, features = head(rownames(enriched.motifs), 6))
    ggsave(paste0(plots.dir, "scatac_motif_fp_", "cluster_", cluster.name, "_", group, ".pdf"), plot = p, width = 8, height = 8)
  }
}

DistanceToTSSInternal <- function(
  peaks,
  genes,
  distance = 200000,
  sep = c("-", "-")
  ) {
  tss <- resize(x = genes, width = 1, fix = 'start')
  genes.extended <- suppressWarnings(
    expr = Extend(
      x = tss, upstream = distance, downstream = distance
    )
  )
  overlaps <- findOverlaps(
    query = peaks,
    subject = genes.extended,
    type = 'any',
    select = 'all'
  )
  hit_matrix <- sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes.extended))
  )
  rownames(x = hit_matrix) <- GRangesToString(grange = peaks, sep = sep)
  colnames(x = hit_matrix) <- genes.extended$gene_name
  return(hit_matrix)
}


CollapseToLongestTranscriptInternal <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}
