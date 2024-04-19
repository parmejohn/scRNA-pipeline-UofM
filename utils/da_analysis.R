#source(paste0(dirname(dirname(dirname(getwd()))),"/utils/misc.R"))
set.seed(333)

CellProportionByCluster <- function(se.integrated) {
  conditions <- unique(se.integrated@meta.data[["group"]])
  group.num <- data.frame(
    . = character(),
    Freq = character(),
    group = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in conditions) {
    cond.num <-
      subset(se.integrated, group == i) %>% Idents() %>% table() %>%
      as.data.frame()
    cond.num$group <- i
    group.num <- rbind(group.num, cond.num)
    # treat.num <- subset(se.integrated, group == ident.2) %>% Idents() %>% table() %>%
    #   as.data.frame()
    # treat.num$group <- ident.2
  }
  
  # group.num <- rbind(ctrl.num, treat.num) %>%
  group.num <- dplyr::rename(group.num, cluster = '.')
  
  cond.prop.plot <- group.num %>%
    group_by(cluster) %>%
    mutate(prop = Freq / sum(Freq)) %>%
    ggplot(aes(x = cluster, y = prop)) +
    geom_col(aes(color = group, fill = group),
             position = position_dodge(0.8),
             width = 0.7) +
    theme(axis.text.x = element_text(angle = 90)) +
    ylim(0, 1)
  PrintSave(cond.prop.plot, 'condition_proportion_per_cluster.pdf')
}

#se.integrated <- readRDS("/home/projects/sc_pipelines/scrna_deanne_harmony_low_res/pipeline/analysis/data/se_integrated_auto_label.rds")
# k = 16
# d = 50
# sample = "sample"
# condition = "group"
# reduced.dims = "INTEGRATED.CCA"
# prop = 0.05
# 
# sc.integrated.milo.traj <- readRDS("/home/projects/sc_pipelines/scrna_deanne_harmony_low_res/pipeline/analysis/data/sc_integrated_milo_traj.rds")

DifferentialAbundanceMilo <-
  function(se.integrated,
           sample,
           condition,
           k,
           d,
           reduced.dims,
           prop = 0.05) {
    DefaultAssay(se.integrated) <- "RNA"
    se.integrated$da.clusters <- Idents(se.integrated)
    #se.integrated@reductions[["harmony"]]@stdev <- as.numeric(apply(se.integrated@reductions[["harmony"]]@cell.embeddings, 2, stats::sd))
    
    # make sure that reduced dims is carried over to save on time
    sc.integrated <-
      as.SingleCellExperiment(se.integrated,  assay = "RNA") # need to convert Seurat to SCE object
    sc.integrated.milo <- Milo(sc.integrated) # make a Milo obj
    
    # Constructing graph -> computes a k-nearest neighbour graph; Graph vertex = single cell, Edges = neighbors
    # built from PCA -> so if using integrated dataset the PCA will be corrected for using MNN
    # sc.integrated.milo.traj <-
    #   buildGraph(
    #     sc.integrated.milo,
    #     k = k,
    #     d = d,
    #     reduced.dim = reduced.dims
    #   ) # d = use the corrected reduction; k is increased or decreased according to plotNhoodSizeHist but usually the same as building knn
    # 
    # # Defining neigborhoods -> group of cells connected by an edge in KNN graph to an index cell
    # sc.integrated.milo.traj <-
    #   makeNhoods(
    #     sc.integrated.milo.traj,
    #     prop = prop,
    #     k = k,
    #     d = d,
    #     reduced_dims = reduced.dims,
    #     refined = TRUE
    #   )
    # # If distribution peak is not ideal, have to change previous k values
    # # can calc bins automatically
    # # how many cells form each neighborhood
    # plotNhoodSizeHist(sc.integrated.milo.traj) # official vignette says it wants the distribution to peak between 50-100 but others say 5 x N samples. Using official vignette's suggestion for now
    # 
    # # counting the cells -> variation in cell numbers between replicates for the same cond to test for DA
    # sc.integrated.milo.traj <-
    #   countCells(sc.integrated.milo.traj,
    #              meta.data = as.data.frame(colData(sc.integrated.milo.traj)),
    #              sample = sample)
    # 
    # # Setting up the design
    # sc.integrated.milo.traj.design <-
    #   data.frame(colData(sc.integrated.milo.traj))[, c(sample, condition)]
    # sc.integrated.milo.traj.design <-
    #   distinct(sc.integrated.milo.traj.design)
    # rownames(sc.integrated.milo.traj.design) <-
    #   sc.integrated.milo.traj.design$sample
    # ## Reorder rownames to match columns of nhoodCounts(milo)
    # sc.integrated.milo.traj.design <-
    #   sc.integrated.milo.traj.design[colnames(nhoodCounts(sc.integrated.milo.traj)), , drop =
    #                                    FALSE]
    # 
    # # store distances -> store distances btw nearest neighbors -> Milo uses this downstream for FDR correction
    # # this step takes long to compute
    # print("starting long compute time")
    # sc.integrated.milo.traj <-
    #   calcNhoodDistance(sc.integrated.milo.traj,
    #                     d = d,
    #                     reduced.dim = reduced.dims)
    # print("end of long compute time")
    # 
    # da.results <-
    #   testNhoods(
    #     sc.integrated.milo.traj,
    #     design = ~ group,
    #     design.df = sc.integrated.milo.traj.design,
    #     reduced.dim = reduced.dims
    #   )
    # 
    ##### testing different method
    ## Build KNN graph neighbourhoods
    milo.obj <- buildGraph(sc.integrated.milo, k=k, d=d, reduced.dim = reduced.dims)
    milo.obj <- makeNhoods(milo.obj, k=k, d=d, refined=T, prop=0.2, refinement_scheme="graph", reduced_dims = reduced.dims)

    ## Count cells in nhoods
    milo.obj <- countCells(milo.obj, samples=sample, meta.data=as.data.frame(colData(milo.obj)))

    # # Setting up the design
    sc.integrated.milo.traj.design <-
      data.frame(colData(milo.obj))[, c(sample, condition)]
    sc.integrated.milo.traj.design <-
      distinct(sc.integrated.milo.traj.design)
    rownames(sc.integrated.milo.traj.design) <-
      sc.integrated.milo.traj.design$sample
    ## Reorder rownames to match columns of nhoodCounts(milo)
    sc.integrated.milo.traj.design <-
      sc.integrated.milo.traj.design[colnames(nhoodCounts(milo.obj)), , drop =
                                       FALSE]
    
    ## Test for differential abundance
    milo_res <- testNhoods(milo.obj, design=~group, design.df=sc.integrated.milo.traj.design, fdr.weighting="graph-overlap", reduced.dim = reduced.dims)
    da.results <- milo_res
    write.table(da.results, "da_diff_test.txt", quote = FALSE,row.names = T, sep = "\t", col.names = T)
    
    sc.integrated.milo.traj <- milo.obj
    
    p1 <-
      ggplot(da.results, aes(PValue)) + 
      geom_histogram(bins = 50) + # should be anti-conservative distribution
      ggtitle("DA Test: Distribution of p-values")
    p2 <-
      ggplot(da.results, aes(logFC,-log10(SpatialFDR))) + #volcano plot to see how many neighborhoods are above a fdr threshold
      geom_point() +
      geom_hline(yintercept = 1) +## Mark significance threshold (10% FDR)
      ggtitle("Significant Neighborhoods over Conditions")
    PrintSave(p1, "milo_pval_distribution.pdf")
    PrintSave(p2, "milo_volcano_plot.pdf")
    
    sc.integrated.milo.traj <-
      buildNhoodGraph(sc.integrated.milo.traj)
    
    # finding DEGs in DA neighborhoods
    da.results <-
      annotateNhoods(sc.integrated.milo.traj, da.results, coldata_col = "da.clusters")
    write.table(da.results, "da_diff_test.txt", quote = FALSE,row.names = T, sep = "\t", col.names = T)
    logcounts(sc.integrated.milo.traj) <-
      log1p(counts(sc.integrated.milo.traj))
    print("buzzbuzz")
    #da.results$da.clusters <- as.numeric(da.results$da.clusters) # if default seurat_clusters?
    
    p4 <- plotDAbeeswarm_fixed(da.results, group.by = "da.clusters") # gives a distribution view instead of UMAP
    if (is.ggplot(p4)) {
    	p4 <- p4 + ggtitle("DA FC Distribution per Cluster")
          PrintSave(p4, "milo_DA_fc_distribution.pdf")
    }
    
    # print("Finding DEGs for DA neighborhoods, this may take a while")
    for (i in levels(droplevels(se.integrated@meta.data[["da.clusters"]]))) {
      print(paste0("Working on DA DE heatmap for cluster ", i))
      logcounts(sc.integrated.milo.traj) <-
        log1p(counts(sc.integrated.milo.traj))
      
      n.da <- sum(da.results$SpatialFDR < 0.1)
      dge_smp <- 0
      if (!is.na(n.da) & n.da == 0) {
        dge_smp <- NULL
      } else {
        # dge_smp <- findNhoodMarkers( ### METHOD IS BEING DEPRECATED
        #   sc.integrated.milo.traj,
        #   da.results,
        #   assay = "logcounts",
        #   gene.offset = FALSE,
        #   da.fdr = 0.1,
        #   aggregate.samples = TRUE,
        #   sample_col = sample,
        #   subset.nhoods = da.results$da.clusters %in% c(i)
        # ) # seeing if this paste method works
        dge_smp <- groupNhoods(
          sc.integrated.milo.traj,
          da.results,
          da.fdr = 0.1,
          subset.nhoods = da.results$da.clusters %in% c(i)
        )
        dge_smp <- findNhoodGroupMarkers(
          sc.integrated.milo.traj,
          dge_smp,
          assay = "logcounts",
          gene.offset = FALSE,
          aggregate.samples = TRUE,
          sample_col = sample,
          subset.nhoods = da.results$da.clusters %in% c(i)
        )
      }
      
      if (!is.null(dge_smp)) {
        #paste("findNhoodMarkers")
        dge.smp.filt <- dge_smp %>%
          filter_at(vars(starts_with("adj.P.Val_")), any_vars(. <= 0.01))
        #paste("filt adjpval")
        print("renaming genes")
        markers <- dge.smp.filt[, "GeneID"]
        #markers <- sub("HLA.", "HLA-", markers) #temporary solution
        print("renamed")
        
        logcounts(sc.integrated.milo.traj) <-
          log1p(counts(sc.integrated.milo.traj))
        print("get log counts again")
        sc.integrated.milo.traj <-
          calcNhoodExpression(sc.integrated.milo.traj, subset.row = markers)
        #paste("traj")
        
        da.results.filt <-
          filter(da.results, da.clusters == i & SpatialFDR < 0.1) # Error in hclust(dist(expr_mat)) : must have n >= 2 objects to cluster -> filtered sets mustve had nothing
        #print("dafilt")
        
        if (!is.null(dim(da.results.filt)[1]) & !is.null(markers)) {
          if (dim(da.results.filt)[1] >= 2 &
              length(markers) >= 2) {
            # have to check if there are any that meet the spatialFDR or marker cutoff to begin with
            print("plotting deg")
            p5 <-
              plotNhoodExpressionDA_fixed(
                sc.integrated.milo.traj,
                da.results.filt,
                features = markers,
                subset.nhoods = da.results$da.clusters %in% c(i),
                assay = "logcounts",
                scale_to_1 = TRUE,
                cluster_features = TRUE,
                show_rownames = FALSE,
                group = c(i)
              ) + 
              ggtitle(paste0(i, ": DA DEGs"))
            PrintSave(p5, paste0("milo_DA_DE_heatmap_", i, ".pdf"))
          }
        }
      }
    }
    # plots the Differential abundance -> uses UMAP reduction from seurat object
    # - nodes = neighborhoods
    #   - color = fold change
    #   - node layout -> based on index node in UMAP of single cells
    # - edges = number cells shared btw adjacent neighborhoods
    p3 <-
      plotNhoodGraphDA(sc.integrated.milo.traj,
                       da.results,
                       layout = "UMAP",
                       alpha = 0.1) +  # for my UMAP negative FC denoted CTRL and positive values denoted TREAT -> will have to look at your original UMAP for density of cells
    # also can look at the your design -> should be ordered correspondingly; https://github.com/MarioniLab/miloR/issues/81
      ggtitle("DA Analysis UMAP")
    
    
    da.results <-
      annotateNhoods(sc.integrated.milo.traj, da.results, coldata_col = "da.clusters")
    #p4 <- plotDAbeeswarm(da.results, group.by = "da.clusters") # gives a distribution view instead of UMAP
    
    #dir.create(paste(plot.path, "da", sep=""))
    #PrintSave(p1, "milo_pval_distribution.pdf")
    #PrintSave(p2, "milo_volcano_plot.pdf")
    PrintSave(p3, "milo_DA_umap.pdf")
    #PrintSave(p4, "milo_DA_fc_distribution.pdf")
    #sc.integrated.milo.traj
    se.integrated$da.clusters <- NULL
    saveRDS(sc.integrated.milo.traj, "sc_integrated_milo_traj.rds")
  }

# x=sc.integrated.milo.traj
# da.res= da.results.filt
# features = markers
# subset.nhoods = da.results$da.clusters %in% c(i)
# assay = "logcounts"
# scale_to_1 = TRUE
# cluster_features = TRUE
# show_rownames = FALSE

plotNhoodExpressionDA_fixed <-
  function(x,
           da.res,
           features,
           alpha = 0.1,
           subset.nhoods = NULL,
           cluster_features = FALSE,
           assay = "logcounts",
           scale_to_1 = FALSE,
           show_rownames = TRUE,
           highlight_features = NULL,
           group = NULL) {
    if (length(features) <= 0 | is.null(features)) {
      stop("features is empty")
    }
    ## Check if features are in rownames(x)
    if (!all(features %in% rownames(x))) {
      stop("Some features are not in rownames(x)")
    }
    ## Check if nhood expression exists
    if (dim(nhoodExpression(x))[2] == 1) {
      warning("Nothing in nhoodExpression(x): computing for requested features...")
      x <-
        calcNhoodExpression(x, assay = assay, subset.row = features)
    }
    ## Check if all features are in nhoodExpression
    if (!all(features %in% rownames(nhoodExpression(x)))) {
      warning("Not all features in nhoodExpression(x): recomputing for requested features...")
      x <-
        calcNhoodExpression(x, assay = assay, subset.row = features)
    }
    
    expr_mat <- nhoodExpression(x)[features,]
    colnames(expr_mat) <- seq_len(ncol(nhoods(x)))
    
    ## Get nhood expression matrix
    if (!is.null(subset.nhoods)) {
      expr_mat <- expr_mat[, subset.nhoods, drop = FALSE]
    }
    
    if (!isFALSE(scale_to_1)) {
      expr_mat <-
        t(apply(expr_mat, 1, function(X)
          (X - min(X)) / (max(X) - min(X))))
      # force NAs to 0?
      if (sum(is.na(expr_mat)) > 0) {
        warning("NA values found - resetting to 0")
        expr_mat[is.na(expr_mat)] <- 0
      }
    }
    
    print("replacing all '-' with '.' because of dataframe conversion issues")
    rownames(expr_mat) <-
      gsub(pattern = "-",
           replacement = ".",
           rownames(expr_mat)) ## To avoid problems when converting to data.frame
    rownames(expr_mat) <-
      gsub(pattern = "\\(",
           replacement = ".",
           rownames(expr_mat)) ## To avoid problems when converting to data.frame
    rownames(expr_mat) <-
      gsub(pattern = "\\)",
           replacement = ".",
           rownames(expr_mat)) ## To avoid problems when converting to data.frame

    pl_df <- data.frame(t(expr_mat)) %>%
      rownames_to_column("Nhood") %>%
      mutate(Nhood = as.double(Nhood)) %>%
      inner_join(da.res, by = "Nhood") %>% # changed from a left join since this will save on compute time
      mutate(logFC_rank = percent_rank(logFC))
    
    pl_df_t <- column_to_rownames(pl_df, "Nhood") %>% select(1:(length(pl_df)-10)) %>% t()
    write.table(pl_df_t, paste0("da_", group, "_markers_by_neighborhood.txt"), quote = FALSE,row.names = T, sep = "\t", col.names = T)
    
    ## Top plot: nhoods ranked by DA log FC
    pl_top <- pl_df %>%
      mutate(is_signif = ifelse(SpatialFDR < alpha, paste0("SpatialFDR < ", alpha), NA)) %>%
      ggplot(aes(logFC_rank, logFC)) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_point(size = 0.2, color = "grey") +
      geom_point(data = . %>% filter(!is.na(is_signif)),
                 aes(color = is_signif),
                 size = 1) +
      theme_bw(base_size = 16) +
      ylab("DA logFC") +
      scale_color_manual(values = "red", name = "") +
      scale_x_continuous(expand = c(0.01, 0)) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
      )
    
    ## Bottom plot: gene expression heatmap
    if (isTRUE(cluster_features)) {
      row.order <- hclust(dist(expr_mat))$order # clustering
      ordered_features <- rownames(expr_mat)[row.order]
    } else {
      ordered_features <- rownames(expr_mat)
    }
    
    # this code assumes that colnames do not begin with numeric values
    # add 'X' to feature names with numeric first characters
    rownames(expr_mat) <-
      str_replace(rownames(expr_mat),
                  pattern = "(^[0-9]+)",
                  replacement = "X\\1")
    
    print("start of error")
    pl_df <- pl_df %>%
      pivot_longer(cols = rownames(expr_mat),
                   names_to = 'feature',
                   values_to = "avg_expr") %>%
      mutate(feature = factor(feature, levels = ordered_features))
    print("past error!")
    
    if (!is.null(highlight_features)) {
      if (!all(highlight_features %in% pl_df$feature)) {
        missing <-
          highlight_features[which(!highlight_features %in% pl_df$feature)]
        warning(
          'Some elements of highlight_features are not in features and will not be highlighted. \nMissing features: ',
          paste(missing, collapse = ', ')
        )
      }
      pl_df <- pl_df %>%
        mutate(label = ifelse(feature %in% highlight_features, as.character(feature), NA))
    }
    
    pl_bottom <- pl_df %>%
      ggplot(aes(logFC_rank, feature, fill = avg_expr)) +
      geom_tile() +
      scale_fill_viridis_c(option = "magma", name = "Avg.Expr.") +
      xlab("Neighbourhoods") + ylab("Features") +
      scale_x_continuous(expand = c(0.01, 0)) +
      theme_classic(base_size = 16) +
      coord_cartesian(clip = "off") +
      theme(
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing = margin(2, 2, 2, 2, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(10, 10, 10, 10)
      )
    
    if (!is.null(highlight_features)) {
      pl_bottom <- pl_bottom +
        geom_text_repel(
          data = . %>%
            filter(!is.na(label)) %>%
            group_by(label) %>%
            summarise(
              logFC_rank = max(logFC_rank),
              avg_expr = mean(avg_expr),
              feature = first(feature)
            ),
          aes(label = label, x = logFC_rank),
          size = 4,
          xlim = c(max(pl_df$logFC_rank) + 0.01, max(pl_df$logFC_rank) + 0.02),
          min.segment.length = 0,
          max.overlaps = Inf,
          seed = 42
        )
      
    }
    
    if (isFALSE(show_rownames)) {
      pl_bottom <- pl_bottom +
        theme(axis.text.y = element_blank())
    }
    
    ## Assemble plot
    (pl_top / pl_bottom) +
      plot_layout(heights = c(1, 4), guides = "collect") &
      theme(legend.justification = c(0, 1),
            legend.margin = margin(0, 0, 0, 50))
  }

# da.res <- da.results
# group.by <- "da.clusters"
# alpha=0.1
# subset.nhoods=NULL

plotDAbeeswarm_fixed <-
  function(da.res,
           group.by = NULL,
           alpha = 0.1,
           subset.nhoods = NULL) {
    if (!is.null(group.by)) {
      if (!group.by %in% colnames(da.res)) {
        stop(
          group.by,
          " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ",
          group.by,
          ")?"
        )
      }
      if (is.numeric(da.res[, group.by])) {
        # stop(group.by, " is a numeric variable. Please bin to use for grouping.")
      }
      da.res <- mutate(da.res, group_by = da.res[, group.by])
    } else {
      da.res <- mutate(da.res, group_by = "g1")
    }
    
    if (!is.factor(da.res[, "group_by"])) {
      message("Converting group_by to factor...")
      da.res <-
        mutate(da.res, group_by = factor(group_by, levels = unique(group_by))) #converts char into numeric for some reason
      #levels(da.res$group_by) <- as.character(unique(da.res.test$group_by))
      # anno_vec <- factor(anno_vec, levels=unique(anno_vec))
    }
    
    if (!is.null(subset.nhoods)) {
      da.res <- da.res[subset.nhoods, ]
    }
    
    
    da.res.test <- da.res %>%
      mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0))
    
    if (length(unique(da.res.test[, "is_signif"])) == 2) {
      # Get position with ggbeeswarm
      beeswarm_pos <- ggplot_build(
        da.res %>%
          mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
          arrange(group_by) %>%
          ggplot(aes(group_by, logFC)) +
          geom_quasirandom()
      )
      
      pos_x <- beeswarm_pos$data[[1]]$x
      pos_y <- beeswarm_pos$data[[1]]$y
      
      n_groups <- unique(da.res$group_by) %>% length()
      
      da.res %>%
        mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
        mutate(logFC_color = ifelse(is_signif == 1, logFC, NA)) %>%
        arrange(group_by) %>%
        mutate(Nhood = factor(Nhood, levels = unique(Nhood))) %>%
        mutate(pos_x = pos_x, pos_y = pos_y) %>%
        ggplot(aes(pos_x, pos_y, color = logFC_color)) +
        scale_color_gradient2() +
        guides(color = "none") +
        xlab(group.by) + ylab("Log Fold Change") +
        scale_x_continuous(breaks = seq(1, n_groups),
                           labels = setNames(levels(da.res$group_by), seq(1, n_groups))) +
        geom_point() +
        coord_flip() +
        theme_bw(base_size = 22) +
        theme(strip.text.y =  element_text(angle = 0))
    } else {
      print("there are no significant spatial FDR for any neighbourhood")
    }
  }
