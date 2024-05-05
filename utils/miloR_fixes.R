set.seed(333)

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
           group = NULL,
           condition = NULL) {
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
    write.table(pl_df_t, paste0("da_", group, "_markers_by_neighborhood_", condition, ".txt"), quote = FALSE,row.names = T, sep = "\t", col.names = T)
    
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