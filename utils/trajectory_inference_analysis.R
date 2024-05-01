set.seed(333)

TrajectoryInferenceSlingshot <- function(se.integrated, start.clus=NULL){
  se.integrated$ti.clusters <- Idents(se.integrated)
  pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))
  
  dim.red <- se.integrated@reductions[["umap"]]@cell.embeddings
  clustering <- se.integrated$ti.clusters
  #counts <- as.matrix(se.integrated@assays[["RNA"]]@layers[["counts"]])
  
  if(is.null(start.clus)){
    
    pdf('ti_no_start_not_smooth.pdf', width = 8, height = 6)
    par(mfrow = c(1, 2))
    plot(dim.red[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
    for (i in levels(clustering)) {
      text(mean(dim.red[clustering == i, 1]), mean(dim.red[clustering == i, 2]), labels = i, font = 1)
    }
    title("scRNA-seq UMAP")
    
    lineages <- getLineages(data = dim.red, clusterLabels = clustering)
    plot(dim.red[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
    lines(as.SlingshotDataSet(lineages), lwd = 3, col = "black")
    title("Lineage Structure")
    graphics.off()

  } else {
    TrajectoryInferenceSlingshotCurved(se.integrated, start.clus)
  }
  fn <- "Rplots.pdf"
  if (file.exists(fn)) {
    #Delete file if it exists
    file.remove(fn)
  }
  se.integrated$ti.clusters <- NULL
}

# https://rnabioco.github.io/cellar/previous/2019/docs/5_trajectories.html
TrajectoryInferenceSlingshotCurved <- function(se.integrated, start.clus){
  # normalized values already found in logcounts in the seurat data
  # log_mat <- log1p(GetAssayData(se.integrated, "RNA"))
  # so <- SetAssayData(se.integrated, "data", new.data = log_mat)
  # test <- as.SingleCellExperiment(so)
  se.integrated$ti.clusters <- Idents(se.integrated)
  #rcl.list <- NA
  #out <- NA
  sce <- as.SingleCellExperiment(se.integrated)
  sce <- suppressWarnings(slingshot(
    sce,
    reducedDim = 'UMAP',
    clusterLabels = 'ti.clusters',
    start.clus = start.clus
  ))
  
  pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))
  dim.red <- se.integrated@reductions[["umap"]]@cell.embeddings
  clustering <- se.integrated$ti.clusters

  pdf('ti_start_smooth.pdf', width = 8, height = 6)
  par(mfrow = c(1, 2))
  plot(dim.red[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
  for (i in levels(clustering)) {
    text(mean(dim.red[clustering == i, 1]), mean(dim.red[clustering == i, 2]), labels = i, font = 1)
  }
  title("scRNA-seq UMAP")
  
  # slo <- SlingshotDataSet(sce) look at the different lineages
  plot(reducedDims(sce)$UMAP, col = brewer.pal(9,'Set1')[sce$ti.clusters], cex = 0.5, pch=16)
  #lines(SlingshotDataSet(sce), lwd=2, col='black')
  line.list <- apply( expand.grid( 1:6, 1:6), 1, paste0, collapse="")
  lineages.list <- c()
  for (i in 1:length(SlingshotDataSet(sce)@curves)){
    lines(SlingshotDataSet(sce)@curves[[paste0("Lineage", i)]], lwd=2, col='black', lty = as.numeric(line.list[i]))
    lineages.list <- append(lineages.list, paste0("Lineage", i))
  }
  legend("bottomright",
         legend=lineages.list,
         col="black",
         lty= as.numeric(line.list[1:3]),
         cex=0.7
         )
  title("Lineage Path Predictions")
  graphics.off()
  
  
  genes_to_test <- VariableFeatures(se.integrated)[1:1000]
  #cnts <- logcounts(sce)[genes_to_test, lineage_cells] # dont need log counts pretty sure
  
  sce.filtered <- sce[genes_to_test,]
  
  # long compute time, this is why the variable features are filtered for
  # can take into account non-normal noise distributions and a greater diversity of non-linear trends
  # takes the counts info with the pseudotime values calculated from slingshot
  print("long compute time for GAM")
  # gam.pval <- apply(cnts, 1, function(z){
  #   d <- data.frame(z = z, 
  #                   ptime = ptime)
  #   tmp <- suppressWarnings(gam(z ~ lo(ptime), data=d))
  #   p <- summary(tmp)[4][[1]][1, 5]
  #   p
  # })
  sce <- fitGAM(sce.filtered) # model the (potentially nonlinear) relationshipships between gene expression and pseudotime
  print("finished GAM computation")
  
  saveRDS(sce, "sce_slingshot.rds")
  # res <- tibble(
  #   id = names(gam.pval),
  #   pvals = gam.pval,
  #   qval = p.adjust(gam.pval, method = "fdr")) %>% 
  #   arrange(qval) # sorts the qval from least to highest
  # rm(gam.pval)
  # gc()
  res <- associationTest(sce) #find if the gene expression is actually associated with pseudotime
  res <- res[complete.cases(res), ]
  res$qval <- p.adjust(res$pvalue, method = "fdr") # BH correction
  
  # iterate over the number of lineages -> output DEG heatmap for each
  for (i in 1:ncol(sce@colData@listData[["slingshot"]])){
    ptime.str <- paste0("slingPseudotime_", i)
    print(paste0("plotting DEGs for TI for ", ptime.str))
    ptime <- sce@colData@listData[[ptime.str]] # pseudotime values
    lineage_cells <- colnames(sce)[!is.na(ptime)] # cells associated with the pseudotime
    ptime <- ptime[!is.na(ptime)]

    
    # get log normalized counts
    to_plot <- NA
    if (length(rownames(res)) >= 100){
      to_plot <- as.matrix(logcounts(sce)[rownames(res[order(res$qval), ])[1:100], lineage_cells]) # get the top 100 genes and filter by
      #to_plot <- as.matrix(logcounts(sce)[res$id, lineage_cells]) 
    } else {
      to_plot <- as.matrix(logcounts(sce)[rownames(res[order(res$qval), ])[1:length(rownames(res))], lineage_cells]) # get the top 100 genes
    }
    
    # arrange cells by pseudotime
    ptime_order <- colnames(to_plot)[order(ptime)]
    
    # add useful annotations
    annotations <- colData(sce)[lineage_cells, 
                                c(ptime.str, 
                                  "ti.clusters")] %>% as.data.frame()
    
    print("start heatmap plotting")
    ha <- HeatmapAnnotation(df = annotations)
    p2 <- Heatmap(to_plot,
                  name = "Log-normalized Counts", 
                  column_title = paste0("DEGs over ", ptime.str),
                  column_order = ptime_order,
                  show_column_names = FALSE,
                  show_row_names = TRUE,
                  top_annotation = ha,
                  row_names_gp = gpar(fontsize = 6),
                  height = 8,
                  width = 8)
    PrintSave(p2, paste0('ti_de_', ptime.str, ".pdf"))
    graphics.off()
    
    ## old method of seperating into clusters, but dont really need since its the same 100 genes
    # p2 <- draw(p2)
    # rcl.list <- row_order(p2)
    # 
    # print("Printing gene + cluster table for DEGs in TI")
    # clu_df <- lapply(names(rcl.list), function(j){
    #   out <- data.frame(GeneID = rownames(as.data.frame(to_plot)[rcl.list[[j]],]), # for some reason rownames cant return a single row value in a matrix; have to convert matrix into a df first
    #                     Cluster = paste0("cluster", j),
    #                     stringsAsFactors = FALSE)
    #   return(out)
    # })  %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
    #   do.call(rbind, .)
    # write.table(clu_df, file= paste0("ti_gene_clusters_", ptime.str, ".txt"), sep="\t", quote=F, row.names=FALSE)
  }
  top.100.gene.list <- NA
  if (length(rownames(res)) >= 100){
    top.100.gene.list <- as.data.frame(rownames(res[order(res$qval), ])[1:100])
    colnames(top.100.gene.list)[1] <- "genes"
  } else {
    top.100.gene.list <- as.data.frame(rownames(res[order(res$qval), ])[1:length(rownames(res))])
    colnames(top.100.gene.list)[1] <- "genes"
  }
  write.table(top.100.gene.list, "top_100_temporally_dynamic_genes.txt", sep="\t", quote=F, row.names=FALSE)
  se.integrated$ti.clusters <- NULL
}
