set.seed(333)

# x=sc.integrated.milo.traj
# da.res= da.results.filt
# features = markers
# subset.nhoods = da.results$da.clusters %in% c(i)
# assay = "logcounts"
# scale_to_1 = TRUE
# cluster_features = TRUE
# show_rownames = FALSE
# alpha = 0.01

plotNhoodExpressionDA_fixed <-
  function(x,
           da.res,
           features,
           alpha = 0.05,
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
        rownames(expr_mat) <-
		gsub(pattern = "\\/",
		replacement = ".",
		rownames(expr_mat)) ## To avoid problems when converting to data.frame
    
    pl_df <- data.frame(t(expr_mat)) %>%
      rownames_to_column("Nhood") %>%
      mutate(Nhood = as.double(Nhood)) %>%
      inner_join(da.res, by = "Nhood") %>% # changed from a left join since this will save on compute time
      mutate(logFC_rank = percent_rank(logFC))
    pl_df_full <- pl_df
    #pl_df_sig <- subset(pl_df, SpatialFDR < alpha)
    
    pl_df_t <- column_to_rownames(pl_df, "Nhood") %>% select(1:(length(pl_df)-10)) %>% t()
    #write.table(pl_df_t, paste0("da_", group, "_markers_by_neighborhood_", condition, ".txt"), quote = FALSE,row.names = T, sep = "\t", col.names = T)
    
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
      ) + 
	ggtitle(paste0(group, ": DA DEGs for ", condition))
    
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
    
#    if (!is.null(highlight_features)) {
#      if (!all(highlight_features %in% pl_df$feature)) {
#        missing <-
#          highlight_features[which(!highlight_features %in% pl_df$feature)]
#        warning(
#          'Some elements of highlight_features are not in features and will not be highlighted. \nMissing features: ',
#          paste(missing, collapse = ', ')
#        )
#      }
#      pl_df <- pl_df %>%
#        mutate(label = ifelse(feature %in% highlight_features, as.character(feature), NA))
#    }
    
#    pl_bottom_old <- pl_df %>%
#      ggplot(aes(logFC_rank, feature, fill = avg_expr)) +
#      geom_tile() +
#      scale_fill_viridis_c(option = "magma", name = "Avg.Expr.") +
#      xlab("Neighbourhoods") + ylab("Features") +
#      scale_x_continuous(expand = c(0.01, 0)) +
#      theme_classic(base_size = 16) +
#      coord_cartesian(clip = "off") +
#      theme(
#        axis.text.x = element_blank(),
#        axis.line.x = element_blank(),
#        axis.ticks.x = element_blank(),
#        axis.line.y = element_blank(),
#        axis.ticks.y = element_blank(),
#        panel.spacing = margin(2, 2, 2, 2, "cm"),
#        legend.margin = margin(0, 0, 0, 0),
#        legend.box.margin = margin(10, 10, 10, 10)
#      )
    
    col_order <- order(pl_df_full$logFC_rank)
    pl_df_sorted <- pl_df_full[col_order,]
    col_order <- unique(pl_df_full[col_order,]$Nhood)
    gap_col_row <- pl_df_sorted[which.min(abs(pl_df_sorted[,"logFC"])),]
    gap_col <- pl_df_sorted[which.min(abs(pl_df_sorted[,"logFC"])),]$Nhood
    gap_col <- match(gap_col, col_order)
    if(gap_col_row$logFC > 0){
                gap_col <- gap_col - 1
    }

    expr_mat_filt <- as.data.frame(expr_mat) %>% select(as.character(unique(da.res$Nhood)))
    slices = diff(c(0, gap_col, ncol(as.matrix(expr_mat_filt))))

    print("start making heatmap")
	is.clustered <- NA
    if (nrow(expr_mat_filt) > 10){
    	pl_bottom <- Heatmap(as.matrix(expr_mat_filt),
        	          name = "Scaled Gene Expression", 
                	  show_column_names = FALSE,
	                  show_row_names = show_rownames,
        	          row_names_gp = gpar(fontsize = 6),
                	  height = 6,
	                  width = 6,
        	          show_column_dend = FALSE,
                	  show_row_dend = FALSE,
	                  column_title = "Neighbourhoods",
        	          column_title_side = "bottom",
                	  km = 10,
	                  show_parent_dend_line = FALSE,
        	          column_order = as.character(col_order),
			  column_split = rep(seq_along(slices), times = slices))
    is.clustered <- TRUE
    } else {
    pl_bottom <- Heatmap(as.matrix(expr_mat_filt),
	                               name = "Scaled Gene Expression",
		                                show_column_names = FALSE,
		                          show_row_names = TRUE,
		                             row_names_gp = gpar(fontsize = 6),
		                          height = 6,
		                          width = 6,
		                          show_column_dend = FALSE,
		                          show_row_dend = FALSE,
		                          column_title = "Neighbourhoods",
		                          column_title_side = "bottom",
		                          column_order = as.character(col_order),
					  row_names_rot = 90,
					  row_names_side = "left",
					  column_split = rep(seq_along(slices), times = slices))			  
    is.clustered <- FALSE
    }
    print("made heatmap")
#    if (!is.null(highlight_features)) {
#      pl_bottom <- pl_bottom +
#        geom_text_repel(
#          data = . %>%
#            filter(!is.na(label)) %>%
#            group_by(label) %>%
#            summarise(
#              logFC_rank = max(logFC_rank),
#              avg_expr = mean(avg_expr),
#              feature = first(feature)
#            ),
#          aes(label = label, x = logFC_rank),
#          size = 4,
#          xlim = c(max(pl_df$logFC_rank) + 0.01, max(pl_df$logFC_rank) + 0.02),
#          min.segment.length = 0,
#          max.overlaps = Inf,
#          seed = 42
#        )
#      
#    }
    
    # if (isFALSE(show_rownames)) {
    #   pl_bottom <- pl_bottom +
    #     theme(axis.text.y = element_blank())
    # }
    
    ## Assemble plot
    # (pl_top / pl_bottom) +
    #   plot_layout(heights = c(1, 4), guides = "collect") &
    #   theme(legend.justification = c(0, 1),
    #         legend.margin = margin(0, 0, 0, 50)
    print("plot_both")
    #print(class(pl_bottom))

    #pl_bottom <- draw(pl_bottom, padding = unit(c(1, 1.25, 1, 1.8), "cm"))
    print("grab expr")
    grob = grid.grabExpr(draw(pl_bottom, padding = unit(c(1, 1.3, 1, 1.8), "cm")))
    print("arrange")

    graphics.off()
    pdf(paste0("milo_DA_DE_heatmap_", group, "_", condition, ".pdf"), width = 8, height = 6)
    pl_both <- grid.arrange(pl_top, grob, nrow = 2, heights = c(1,4))
    graphics.off()

    if(is.clustered){
    print("get row order")
    pl_bottom <- draw(pl_bottom)
    rcl.list <- row_order(pl_bottom)

    print("Printing gene + cluster table for DEGs in miloR")
    clu_df <- lapply(names(rcl.list), function(p){
      out <- data.frame(GeneID = rownames(as.data.frame(expr_mat_filt)[rcl.list[[p]],]), # for some reason rownames cant return a single row value in a matrix; have to convert matrix into a df first
                        Cluster = paste0("cluster", p),
                        stringsAsFactors = FALSE)
      return(out)
    })  %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
      do.call(rbind, .)
    write.table(clu_df, file = paste0("da_", group, "_markers_by_neighborhood_", condition, ".txt"), sep="\t", quote=F, row.names=FALSE)
    }
    print("writing table")
    write.table(expr_mat_filt, paste0("da_", group, "_markers_by_neighborhood_", condition,"_expr_matrix.txt"), quote = FALSE,row.names = T, sep = "\t", col.names = T)

    return(pl_both)
  }


# da.res <- da.results
# group.by <- "da.clusters"
# alpha=0.1
# subset.nhoods=NULL

plotDAbeeswarm_fixed <-
  function(da.res,
           group.by = NULL,
           alpha = 0.05,
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

calcNhoodExpression <- function(x, assay="logcounts", subset.row=NULL, exprs=NULL){

	    if(is(x, "Milo")){
		            # are neighbourhoods calculated?
		            if(ncol(nhoods(x)) == 1 & nrow(nhoods(x)) == 1){
				                stop("No neighbourhoods found - run makeNhoods first")
        }

        if(!is.null(assay(x, assay))){
		            n.exprs <- .calc_expression_fixed(nhoods=nhoods(x),
							                                        data.set=assay(x, assay),
												                                        subset.row=subset.row)
	            nhoodExpression(x) <- n.exprs
		                return(x)
		            }
	    } else if(is(x, "Matrix")){
		            if(is.null(exprs)){
				                stop("No expression data found. Please specific a gene expression matrix to exprs")
	            } else{
			                n.exprs <- .calc_expression(nhoods=x,
								                                            data.set=exprs,
													                                            subset.row=subset.row)
		                x.milo <- Milo(SingleCellExperiment(assay=list(logcounts=exprs)))
				            nhoodExpression(x.milo) <- n.exprs
				            return(x.milo)
					            }
	        }
}

.calc_expression_fixed <- function(nhoods, data.set, subset.row=NULL){
  # neighbour.model <- matrix(0L, ncol=length(nhoods), nrow=ncol(data.set))
  #
  # for(x in seq_along(1:length(nhoods))){
  #     neighbour.model[nhoods[[x]], x] <- 1
  # }
  
  if(!is.null(subset.row)){
    if(is(data.set[subset.row,], "Matrix")){
	    print(dim(Matrix::t(nhoods)))
	    print(dim(data.set[subset.row,]))
      neigh.exprs <- Matrix::tcrossprod(Matrix::t(nhoods), data.set[subset.row,])
    }else{
	print(dim(Matrix::t(nhoods)))
    	print(dim(as(as(as(data.set[subset.row,], "dMatrix"), "generalMatrix"), "unpackedMatrix")))
	#neigh.exprs <- Matrix::tcrossprod(Matrix::t(nhoods), as(data.set[subset.row,], "sparseMatrix"))
      neigh.exprs <- Matrix::tcrossprod(Matrix::t(nhoods), t(as(as(as(data.set[subset.row,], "dMatrix"), "generalMatrix"), "unpackedMatrix")))
      #neigh.exprs <- Matrix::tcrossprod(Matrix::t(nhoods), as(data.set[subset.row,], "dgeMatrix"))
    }
  } else{
    neigh.exprs <- Matrix::tcrossprod(Matrix::t(nhoods), data.set)
  }
  neigh.exprs <- t(apply(neigh.exprs, 2, FUN=function(XP) XP/colSums(nhoods)))
  
  if(is.null(subset.row)){
    rownames(neigh.exprs) <- rownames(data.set)
  } else{
    rownames(neigh.exprs) <- rownames(data.set[subset.row, , drop=FALSE])
  }
  
  return(neigh.exprs)
}

## this function will take into account which genes has a high enough logFC and adj.Pval after performing DEGs on the group of neighborhoods in a given cluster
dynamic_filter_function <- function(data, logfc_pattern = "^logFC_", pval_pattern = "^adj\\.P.Val_", logfc_threshold = 1, pval_threshold = 0.01) {
  logfc_cols <- grep(logfc_pattern, names(data), value = TRUE)
  pval_cols <- grep(pval_pattern, names(data), value = TRUE)
  
  filter_expr <- purrr::map2(logfc_cols, pval_cols, function(logfc_col, pval_col) {
    rlang::expr((abs(!!rlang::sym(logfc_col)) >= logfc_threshold & !!rlang::sym(pval_col) <= pval_threshold))
  }) %>%
    purrr::reduce(~rlang::expr((!!.x) | (!!.y)))
  #print(filter_expr)
  
  data %>%
    filter(!!filter_expr)
}

DAGseaComparison <- function(de.markers, cluster.name, group, fgsea.sets){
  
  cluster.genes <- de.markers %>%
    arrange(desc(avg_logFC)) %>% 
    dplyr::select(gene, avg_logFC) # use avg_log2FC as ranking for now; https://www.biostars.org/p/9526168/
  
  ranks <- deframe(cluster.genes)
  
  fgseaRes <- fgsea(fgsea.sets, stats = ranks)
  fgseaRes$leadingEdge <- as.character(fgseaRes$leadingEdge)
  write.table(fgseaRes, paste0("milo_gsea_cluster_", cluster.name, "_", group,".txt"), quote = FALSE,row.names = T, sep = "\t", col.names = T)
  fgseaRes <- filter(fgseaRes, padj <= 0.05 & size >= 3) %>% arrange(desc(NES)) # be more lenient with the pval cutoff since exploratory analysis
  fgseaRes$Enrichment = ifelse(fgseaRes$NES > 0, "Up-regulated", "Down-regulated") 

  filtRes <-  rbind(head(fgseaRes, n = 10),
                    tail(fgseaRes, n = 10))
  
  if (dim(filtRes)[1] > 0){
    p <- ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
      geom_point( aes(fill = Enrichment, size = size), shape=21) + # size is equal to the amount of genes found in the given gene set
      scale_size_continuous(range = c(2,10)) +
      geom_hline(yintercept = 0) +
      coord_flip() +
      theme_bw() +
      labs(x="Pathway", y="Normalized Enrichment Score") + 
      ggtitle(paste0("GSEA: ", cluster.name, " ", group))
    
    #dir.create(paste(plot.path, 'gsea', sep=''))
    PrintSave(p, paste0("milo_gsea_cluster_", cluster.name, "_", group, '.pdf'), w=12)
  } else {
    print("No pathways are particularly enriched")
  }
}
