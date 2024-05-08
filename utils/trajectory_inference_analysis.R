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
  
  ## If the main condition has more than 2 variables -> perform the conditionTest on it
  ## This will allow one to find temporally DEGs between 2 conditions; have not tested if this method works with 3 variables
  if (length(unique(sce@colData@listData[["group"]])) >=2){
    sce <- fitGAM(sce.filtered, conditions = factor(sce@colData@listData[["group"]])) # model the (potentially nonlinear) relationshipships between gene expression and pseudotime
    
    ## Find temporally different genes between conditions
    ## Unsure if the condition or association test will work if the main condition has > 2 variables to work with -> will have to test
    condRes <- conditionTest(sce, l2fc = log2(2))
    condRes$padj <- p.adjust(condRes$pvalue, "fdr")
    conditionGenes <- rownames(condRes)[condRes$padj <= 0.05]
    conditionGenes <- conditionGenes[!is.na(conditionGenes)]
    
    # provide smoothing to time series -> rmv variation btw time steps -> rmv noise and get the signal of underlying causes
    # did so for the conditionTest to condense the cell and pseudotime information -> would end up with too many columns to visualize
    # scale smoothed time series data https://hectorrdb.github.io/condimentsPaper/articles/TGFB.html
    conditionGenes.smooth <- predictSmooth(sce, gene = conditionGenes, nPoints = 100, tidy = FALSE) %>%
      log1p()
    conditionGenes.smooth.scaled <- t(apply(conditionGenes.smooth,1, scales::rescale))
    
    mainconditions.list <- list()
    c = 1
    for (i in unique(sce@colData@listData[["group"]])){
      
      ## dynamically create pheatmaps and stitch together
      firstcol = grep(paste0("condition", i, "_point1$"), colnames(conditionGenes.smooth))[1]
      lastcol = grep(paste0("condition", i, "_point100$"), 
                     colnames(conditionGenes.smooth))[length(grep(paste0("condition", i, "_point100$"), colnames(conditionGenes.smooth)))]
      
      cond.pheatmap <- NA
      if (length(mainconditions.list) < 1){
        cond.pheatmap <- pheatmap::pheatmap(conditionGenes.smooth.scaled[,  firstcol:lastcol],
                                            cluster_cols = FALSE,
                                            show_rownames = FALSE, show_colnames = FALSE, main = i, legend = FALSE, silent = TRUE)
      } else if (c < length(unique(sce@colData@listData[["group"]]))){
        cond.pheatmap <- pheatmap::pheatmap(conditionGenes.smooth.scaled[mainconditions.list[[1]][["order"]], firstcol:lastcol],
                                            cluster_cols = FALSE, cluster_rows = FALSE,
                                            show_rownames = FALSE, show_colnames = FALSE, main = i, legend = FALSE, silent = TRUE)
      } else if (c == length(unique(sce@colData@listData[["group"]]))){
        cond.pheatmap <- pheatmap::pheatmap(conditionGenes.smooth.scaled[mainconditions.list[[1]][["order"]], firstcol:lastcol],
                                            cluster_cols = FALSE, cluster_rows = FALSE,
                                            show_rownames = FALSE, show_colnames = FALSE, main = i, legend = TRUE, silent = TRUE)
      }
      mainconditions.list <- append(mainconditions.list, cond.pheatmap)
      c = c + 1
    }
    
    hm.coord <- c()
    for (i in 1:length(mainconditions.list)){
      if (i %% 4 != 0){
        hm.coord <- c(hm.coord, i)
      } 
    }
    mainconditions.list.filt <- mainconditions.list[-hm.coord]
    
    pdf(paste0('ti_deg_between_', "group", ".pdf"), width = 8, height = 6)
    do.call("grid.arrange", c(mainconditions.list.filt, ncol=length(mainconditions.list.filt), top="Temporally DEGs across Main Grouping"))
    #PrintSave(p3, paste0('ti_deg_between_', "group", ".pdf"))
    graphics.off()
    
    write.table(as.data.frame(conditionGenes), "ti_de_between_group.txt", sep="\t", quote=F, row.names=FALSE)
    
  } else {
    sce <- fitGAM(sce.filtered)
  }
  print("finished GAM computation")
  
  saveRDS(sce, "sce_slingshot.rds")
  
  # res <- tibble(
  #   id = names(gam.pval),
  #   pvals = gam.pval,
  #   qval = p.adjust(gam.pval, method = "fdr")) %>% 
  #   arrange(qval) # sorts the qval from least to highest
  # rm(gam.pval)
  # gc()
  
  res <- associationTest(sce, l2fc = log2(2), lineages = TRUE) # find if the gene expression is associated with pseudotime for each specific lineage comparison
  #res <- res[complete.cases(res), ]
  #res$qval <- p.adjust(res$pvalue, method = "fdr") # BH correction
  
  # iterate over the number of lineages -> output how the top 100 DEGs are changing over each different lineage
  for (i in 1:ncol(sce@colData@listData[["slingshot"]])){
    ptime.str <- paste0("slingPseudotime_", i)
    print(paste0("plotting DEGs for TI for ", ptime.str))
    ptime <- sce@colData@listData[[ptime.str]] # pseudotime values
    lineage_cells <- colnames(sce)[!is.na(ptime)] # cells associated with the pseudotime
    ptime <- ptime[!is.na(ptime)]
    
    lineage.cols = grep(paste0("lineage", i), colnames(res))
    lineage.res <- res[, c(lineage.cols)]
    lineage.res <- lineage.res[complete.cases(lineage.res), ]
    
    lineage.pval.cols = grep(paste0("pvalue"), colnames(lineage.res))
    
    lineage.qval.colnames <- c()
    for (j in lineage.pval.cols){
      qval.col <- paste0(colnames(lineage.res)[j],"_qval")
      lineage.res[, qval.col] <-  p.adjust(lineage.res[,j], method = "fdr")
      
      lineage.qval.colnames <- c(lineage.qval.colnames, qval.col)
    }
    
    lineage.res <-  lineage.res %>% select(all_of(lineage.qval.colnames))
    lineage.res <- filter_all(lineage.res, any_vars(. < 0.05))
    
    ## Rank genes by lowest FDR value
    lineage.res <- mutate(lineage.res, Rank = min_rank(paste(eval(parse(text = lineage.qval.colnames)))))
    
    # get log normalized counts
    to_plot <- NA
    if (length(rownames(lineage.res)) >= 100){
      to_plot <- as.matrix(logcounts(sce)[rownames(lineage.res[order(lineage.res$Rank), ])[1:100], lineage_cells]) # get the top 100 genes and filter by
      #to_plot <- as.matrix(logcounts(sce)[res$id, lineage_cells]) 
    } else {
      to_plot <- as.matrix(logcounts(sce)[rownames(lineage.res[order(lineage.res$Rank), ])[1:length(rownames(lineage.res))], lineage_cells]) # get the top 100 genes
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
    
    ## old method of seperating into clusters, but dont really need since cluster names are printed
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
    write.table(lineage.res, paste0("ti_DEGs_qval_full_lineage_", i,".txt"), sep="\t", quote=F, row.names=FALSE)
    
  }
  # top.100.gene.list <- NA
  # if (length(rownames(res)) >= 100){
  #   top.100.gene.list <- as.data.frame(rownames(lineage.res[order(lineage.res$qval), ])[1:100])
  #   colnames(top.100.gene.list)[1] <- "genes"
  # } else {
  #   top.100.gene.list <- as.data.frame(rownames(lineage.res[order(lineage.res$qval), ])[1:length(rownames(lineage.res))])
  #   colnames(top.100.gene.list)[1] <- "genes"
  # }

  #write.table(top.100.gene.list, paste0("ti_DEGs_qval_full_lineage_", i,".txt"), sep="\t", quote=F, row.names=FALSE)
  se.integrated$ti.clusters <- NULL
}
