source("./utils/misc.R")

# code source: https://hbctraining.github.io/scRNA-seq_online/lessons/09_merged_SC_marker_identification.html
GetConserved <- function(cluster){ # replace or edit function?, 
  FindConservedMarkers(se.integrated,
                       ident.1 = cluster,
                       grouping.var = "group", # condition should be named group by default in pipeline
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}

IdentifyCellMarkers <- function(se.integrated, outdir, plot.path){
  se.markers.presto <- RunPrestoAll(se.integrated, only.pos = T, logfc.threshold = 0.1, min.pct = 0.01)
  write.table(se.markers.presto, paste(outdir, "/se_markers_presto_integrated.txt", sep = ''), 
              quote = FALSE,row.names = T, sep = "\t", col.names = T)
  
  se.markers.presto.top3 <- as.data.frame(se.markers.presto %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))
  expr.heatmap <- DoHeatmap(se.integrated, features = se.markers.presto.top3$gene)
  PrintSave(expr.heatmap, 'top3_markers_expr_heatmap.pdf', plot.path)
  
  ##### Conserved markers across the conditions #####
  if (length(unique(se.integrated@meta.data[["group"]])) == 2){
    cluster.num <- nlevels(se.integrated@meta.data[["seurat_clusters"]])-1
    conserved.markers <- map_dfr(0:cluster.num, GetConserved) # cluster 27 doesnt show any conserved because its purely made of 1 condtion (CTRL)
    conserved.markers.top.2 <- conserved.markers %>% 
      mutate(avg_fc = (eval(as.symbol(paste(unique(se.integrated@meta.data[["group"]])[1], '_avg_log2FC', sep=''))) + 
                         eval(as.symbol(paste(unique(se.integrated@meta.data[["group"]])[2], '_avg_log2FC', sep='')))) /2) %>% # extract from group name
      group_by(cluster_id) %>% 
      top_n(n = 2, wt = avg_fc)
    
    # color intensity denotes average expression across all cells in a class
    conserved.markers.dotplot <- DotPlot(se.integrated, features = conserved.markers.top.2$gene, cols = c("blue", "red"), dot.scale = 8, split.by = "group") +
      RotatedAxis()
    PrintSave(conserved.markers.dotplot, 'conserved_marker_unlabelled.pdf', plot.path, w=20, h = 16)
  }
}

ReferenceMarkerMapping <- function(reference, query, dims, plot.path){
  reference <- NormalizeData(reference)
  reference <- FindVariableFeatures(object = reference)
  
  query <- NormalizeData(query)
  query <- FindVariableFeatures(query)
  
  #https://satijalab.org/seurat/reference/transferdata
  se.anchors <- FindTransferAnchors(reference = reference, query = query, dims = 1:dims)
  reference$reference.cell.meta <- Idents(reference)
  se.predictions <- TransferData(anchorset = se.anchors, refdata = reference$reference.cell.meta, dims = 1:dims)
  query <- AddMetaData(query, metadata = se.predictions)
  
  prediction.scores <- query@meta.data[, grepl("^prediction|RNA_snn_res.1", names(query@meta.data))]
  prediction.scores <- prediction.scores[,-which(names(prediction.scores) == "prediction.score.max")]
  colnames(prediction.scores) <- gsub("prediction.score.", "", colnames(prediction.scores))
  prediction.scores <- melt(prediction.scores, id.vars = "RNA_snn_res.1", variable.name = "source", value.name = "score")
  
  prediction.matrix <- tapply(prediction.scores$score, list(prediction.scores$RNA_snn_res.1, prediction.scores$source), median)
  se.hm <- pheatmap(prediction.matrix, cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(c("white","red"))(200), display_numbers = FALSE, silent = TRUE)
  PrintSave(se.hm, "reference_marker_mapping_heatmap.pdf", plot.path)
  
  cluster.names <- colnames(prediction.matrix)[max.col(prediction.matrix,ties.method="first")]
  names(cluster.names) <- levels(query)
  query <- RenameIdents(query, cluster.names)
  query$celltype <- Idents(query)
  # query$celltype <- factor(query$celltype, levels = c(cluster.names))
  Idents(query) <- 'celltype'
  query$ident <- query$celltype
  
  umap.labelled <- DimPlot(query, reduction = "umap", group.by = "celltype")
  PrintSave(umap.labelled, "integrated_umap_labelled.pdf", plot.path)
  reference$reference.cell.meta <- NULL
  
  query
}
