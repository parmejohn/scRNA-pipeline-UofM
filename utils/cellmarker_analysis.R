set.seed(333)

# code source: https://hbctraining.github.io/scRNA-seq_online/lessons/09_merged_SC_marker_identification.html
GetConserved <- function(cluster){ # replace or edit function?, 
  FindConservedMarkers(se.integrated,
                       ident.1 = cluster,
                       grouping.var = "group", # condition should be named group by default in pipeline
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

IdentifyCellMarkers <- function(se.integrated, plots.format){
  se.markers.presto <- RunPrestoAll(se.integrated, only.pos = T, logfc.threshold = 0.1, min.pct = 0.01, min.cells.group = 1)
  write.table(se.markers.presto, paste("se_markers_presto_integrated.txt", sep = ''), 
              quote = FALSE,row.names = T, sep = "\t", col.names = T)
  
  se.markers.presto.top3 <- as.data.frame(se.markers.presto %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
  expr.heatmap <- DoHeatmap(subset(se.integrated, downsample = 300), 
                            features = se.markers.presto.top3$gene, 
                            assay = "RNA", slot = "data", label = FALSE) + # will only plot a max of 300 cells per identity; ggplot has a limit of 30,000 cells total
    ggtitle("Top 3 Cell Markers per Cluster")
  PrintSave(expr.heatmap, 'top3_markers_expr_heatmap', plots.format)
  PrintSave(expr.heatmap, 'top3_markers_expr_heatmap', "jpeg")
  
  ##### Conserved markers across the conditions #####
  #   Error in `levels<-`(`*tmp*`, value = as.character(levels)) : -> fix later
  # factor level [32] is duplicated
  # Calls: IdentifyCellMarkers -> DotPlot -> factor
  
  if (length(unique(se.integrated@meta.data[["group"]])) == 2){
    cluster.num <- nlevels(se.integrated@meta.data[["seurat_clusters"]])-1
    conserved.markers <- map_dfr(0:cluster.num, GetConserved) # cluster 27 doesnt show any conserved because its purely made of 1 condtion (CTRL)
    conserved.markers.top.2 <- conserved.markers %>%
      mutate(avg_fc = (eval(as.symbol(paste0(unique(se.integrated@meta.data[["group"]])[1], '_avg_log2FC'))) +
                         eval(as.symbol(paste0(unique(se.integrated@meta.data[["group"]])[2], '_avg_log2FC')))) /2) %>%
      group_by(cluster_id) %>%
      top_n(n = 2, wt = avg_fc)
    write.table(conserved.markers.top.2, paste("se_markers_conserved_integrated.txt", sep = ''), 
                quote = FALSE,row.names = T, sep = "\t", col.names = T)

    # color intensity denotes average expression across all cells in a class
    conserved.markers.dotplot <- DotPlot(se.integrated, features = unique(conserved.markers.top.2$gene), cols = c("blue", "red"), dot.scale = 8, split.by = "group") +
      RotatedAxis() + 
      ggtitle("Conserved Markers across conditions per Cluster") + 
      xlab("Genes") +
      ylab("Cluster and Condition Name")
    PrintSave(conserved.markers.dotplot, 'conserved_marker_unlabelled', width=20, height = 16, plots.format = plots.format)
    PrintSave(conserved.markers.dotplot, 'conserved_marker_unlabelled', width=20, height = 16, plots.format = "jpeg")
  }
}

ReferenceMarkerMapping <- function(reference, query, dims, plots.format){
  DefaultAssay(reference) <- "RNA"
  DefaultAssay(query) <- "RNA"
  reference <- NormalizeData(reference)
  reference <- FindVariableFeatures(object = reference)
  reference <- ScaleData(object = reference)
  
  query <- NormalizeData(query)
  query <- FindVariableFeatures(query)
  
  #https://satijalab.org/seurat/reference/transferdata
  se.anchors <- FindTransferAnchors(reference = reference, query = query, dims = 1:dims)
  reference$reference.cell.meta <- Idents(reference)
  se.predictions <- TransferData(anchorset = se.anchors, refdata = reference$reference.cell.meta, dims = 1:dims)
  query <- AddMetaData(query, metadata = se.predictions)
  
  prediction.scores <- query@meta.data[, grepl("^prediction|snn_res", names(query@meta.data))]
  meta.data.snn <- names(se.integrated@meta.data)[grepl("^prediction|snn_res", names(se.integrated@meta.data))]
  prediction.scores <- prediction.scores[,-which(names(prediction.scores) == "prediction.score.max")]
  colnames(prediction.scores) <- gsub("prediction.score.", "", colnames(prediction.scores))
  prediction.scores <- melt(prediction.scores, id.vars = meta.data.snn, variable.name = "source", value.name = "score")
  
  prediction.matrix <- tapply(prediction.scores$score, list(prediction.scores[,meta.data.snn], prediction.scores$source), median)
  se.hm <- pheatmap(prediction.matrix, cluster_rows = FALSE, cluster_cols = FALSE, 
                    color = colorRampPalette(c("white","red"))(200), display_numbers = FALSE, 
                    main = "Reference Marker Prediction Scores", silent = TRUE)
  PrintSave(se.hm, "reference_marker_mapping_heatmap", plots.format)
  PrintSave(se.hm, "reference_marker_mapping_heatmap", "jpeg")
  
  cluster.names <- colnames(prediction.matrix)[max.col(prediction.matrix,ties.method="first")]
  names(cluster.names) <- levels(query)
  query <- RenameIdents(query, cluster.names)
  query$celltype <- Idents(query)
  # query$celltype <- factor(query$celltype, levels = c(cluster.names))
  Idents(query) <- 'celltype'
  query$ident <- query$celltype
  
  umap.labelled <- DimPlot(query, reduction = "umap", group.by = "celltype", label = TRUE, alpha = 0.5) +
    ggtitle("UMAP Reference Labelled Clusters")
  PrintSave(umap.labelled, "integrated_umap_labelled", plots.format)
  PrintSave(umap.labelled, "integrated_umap_labelled", "jpeg")
  
  p2 <- dittoBarPlot(
    object = query,
    var = "celltype",
    group.by = "group",
    retain.factor.levels=TRUE,
    main = "Percent of cells by Labelled Clusters")
  PrintSave(p2, "percent_cells_group_labelled", plots.format)
  PrintSave(p2, "percent_cells_group_labelled", "jpeg")
  
  reference$reference.cell.meta <- NULL
  
  query
}
