set.seed(333)

TrajectoryInferenceSlingshot <- function(se.integrated, start.clus=NULL){
  se.integrated$ti.clusters <- Idents(se.integrated)
  pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))
  
  dim.red <- se.integrated@reductions[["umap"]]@cell.embeddings
  clustering <- se.integrated$ti.clusters
  counts <- as.matrix(se.integrated@assays[["RNA"]]@layers[["counts"]])
  
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
    # p1 <- recordPlot()
    # plot.new()
    # PrintSave(p1, 'ti_no_start_not_smooth.pdf')
    # if(dev.cur() > 1){
    #   par(mfrow=c(1,1))
    #   dev.off()
    # }

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
TrajectoryInferenceSlingshotCurved <- function(se.integrated, start.clus, km=10){
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
  lines(SlingshotDataSet(sce), lwd=2, col='black')
  title("Lineage Path Predictions")
  graphics.off()
  # p1 <- recordPlot()
  # plot.new()
  # PrintSave(p1,'ti_start_smooth.pdf')
  # if(dev.cur() > 1){
  #   par(mfrow=c(1,1))
  #   dev.off()
  # }
  
  saveRDS(sce, "sce_slingshot.rds")
  
  
  # iterate over the number of lineages -> output DEG heatmap for each
  for (i in 1:ncol(sce@colData@listData[["slingshot"]])){
    ptime.str <- paste0("slingPseudotime_", i)
    print(paste0("plotting DEGs for TI for ", ptime.str))
    ptime <- sce@colData@listData[[ptime.str]]
    lineage_cells <- colnames(sce)[!is.na(ptime)]
    ptime <- ptime[!is.na(ptime)]
    genes_to_test <- VariableFeatures(se.integrated)[1:500]
    cnts <- logcounts(sce)[genes_to_test, lineage_cells] # dont need log counts pretty sure
    
    # long compute time, this is why the variable features are filtered for
    # can take into account non-normal noise distributions and a greater diversity of non-linear trends
    # takes the counts info with the pseudotime values calculated from slingshot
    print("long compute time for GAM")
    gam.pval <- apply(cnts, 1, function(z){
      d <- data.frame(z = z, 
                      ptime = ptime)
      tmp <- suppressWarnings(gam(z ~ lo(ptime), data=d))
      p <- summary(tmp)[4][[1]][1, 5]
      p
    })
    
    res <- tibble(
      id = names(gam.pval),
      pvals = gam.pval,
      qval = p.adjust(gam.pval, method = "fdr")) %>% 
      arrange(qval)
    
    # get log normalized counts 
    #to_plot <- as.matrix(logcounts(sce)[res$id[1:100], lineage_cells]) #limited to the first 100 genes which will miss biological importance
    to_plot <- as.matrix(logcounts(sce)[res$id, lineage_cells]) #limited to the first 100 genes which will miss biological importance
    
    # arrange cells by pseudotime
    ptime_order <- colnames(to_plot)[order(ptime)]
    
    # add useful annotations
    annotations <- colData(sce)[lineage_cells, 
                                c(ptime.str, 
                                  "ti.clusters")] %>% as.data.frame()
    
    ha <- HeatmapAnnotation(df = annotations)
    p2 <- Heatmap(to_plot,
                  name = "Log-normalized Counts", 
                  column_title = paste0("DEGs over ", ptime.str),
                  km = km,
                  column_order = ptime_order,
                  show_column_names = FALSE,
                  show_row_names = FALSE,
                  top_annotation = ha)
    PrintSave(p2, paste0('ti_de_', ptime.str, ".pdf"))
    
    p2 <- draw(p2)
    rcl.list <- row_order(p2)
    
    print("Printing gene + cluster table for DEGs in TI")
    # print(lapply(rcl.list, function(x) length(x)))
    # print(class(to_plot))
    # print(sum(duplicated(rownames(to_plot))))
    
    # for (j in 1:length(rcl.list)){
    #   if (j == 1) {
    #     clu <- t(t(row.names(to_plot[rcl.list[[j]],])))
    #     out <- cbind(clu, paste("cluster", j, sep=""))
    #     colnames(out) <- c("GeneID", "Cluster")
    #   } else {
    #     clu <- t(t(row.names(to_plot[rcl.list[[j]],])))
    #     clu <- cbind(clu, paste("cluster", j, sep=""))
    #     out <- rbind(out, clu)
    #   }
    # }
    
    clu_df <- lapply(names(rcl.list), function(j){
      out <- data.frame(GeneID = rownames(as.data.frame(to_plot)[rcl.list[[j]],]), # for some reason rownames cant return a single row value in a matrix; have to convert matrix into a df first
                        Cluster = paste0("cluster", j),
                        stringsAsFactors = FALSE)
      return(out)
    })  %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
      do.call(rbind, .)
    write.table(clu_df, file= paste0("ti_gene_clusters_", ptime.str, ".txt"), sep="\t", quote=F, row.names=FALSE)
  }
  se.integrated$ti.clusters <- NULL
}
