set.seed(333)

PreprocessingSeurat <- function(seurat_object){
  DefaultAssay(seurat_object) <- "RNA"
  
  seurat_object <- NormalizeData(seurat_object) # feature counts for each cell is divided by total counts for the cell; not based on sample
  seurat_object <- FindVariableFeatures(seurat_object) # chooses top variable features (genes); all cells are merged together anyways
  
  # regressing percent.mt for now, for cell cycle could hold importance for cells undergoing differentiation
  # also automating whether the difference in clustering is affected by cell cycle phase requires visual checks
  seurat_object <- ScaleData(seurat_object, vars.to.regress = c("percent.mt")) # scale and center the data, while regressing specified features
  seurat_object <- RunPCA(seurat_object) # dimensional reduction -> large set of genes into smaller metafeature clusters
  
  if ("ATAC" %in% SeuratObject::Assays(seurat_object)) {
    DefaultAssay(seurat_object) <- "ATAC"
    seurat_object <- RunTFIDF(seurat_object, assay = "ATAC")
    print("runtfidf")
    
    seurat_object <- FindTopFeatures(seurat_object, min.cutoff = "q5", assay = "ATAC")
    print("found top features")
    
    seurat_object <- RunSVD(seurat_object, assay = "ATAC")
    print("run svd")
    
    DefaultAssay(seurat_object) <- "RNA"
    print("reset default")
    print(seurat_object)
  }
  return(seurat_object)
}

# integratation
IntegrateSamples <- function(seurat_obj_list, reduction){
  if (length(seurat_obj_list) >= 2){
    se.merged <- merge(seurat_obj_list[[1]], seurat_obj_list[c(2:length(seurat_obj_list))])
    se.merged[["RNA"]] <- split(se.merged[["RNA"]], f = se.merged$group) # whether the data is split by sample or by treatment does not matter for downstream analysis

    se.merged <- UpdateSeuratObject(se.merged)
    se.merged.preprocessed <- PreprocessingSeurat(se.merged)
    
    if ("ATAC" %in% SeuratObject::Assays(se.merged.preprocessed)) {
      se.merged.preprocessed.atac <- IntegrateAtac(se.merged.preprocessed)
      se.merged.preprocessed[['ATAC']] <- se.merged.preprocessed.atac[['ATAC']]
      se.merged.preprocessed[['integrated.lsi']] <- se.merged.preprocessed.atac[['integrated.lsi']]
    }
    
    # CCA integration background; https://hbctraining.github.io/scRNA-seq_online/lessons/06_integration.html
    # - It is a form of PCA, in that it identifies the greatest sources of variation in the data, but only if it is shared or conserved across the conditions/groups
    DefaultAssay(se.merged.preprocessed) <- "RNA"
    print(se.merged.preprocessed)
    
    if (reduction == "harmony"){
      se.integrated <- IntegrateLayers(object = se.merged.preprocessed, method = HarmonyIntegration,
                                       orig.reduction = "pca", new.reduction = "harmony",
                                       assay = "RNA",
                                       verbose = T)
    } else if (reduction == "integrated.cca"){
      se.integrated <- IntegrateLayers(object = se.merged.preprocessed, method = CCAIntegration,
                                       orig.reduction = "pca", new.reduction = "integrated.cca",
                                       assay = "RNA",
                                       verbose = T)
    } else if (reduction == "integrated.mnn"){
      se.integrated <- IntegrateLayers(object = se.merged.preprocessed, method = FastMNNIntegration,
                                       orig.reduction = "pca", new.reduction = "integrated.mnn",
                                       assay = "RNA",
                                       verbose = T)
    } else {
      stop(print0(reduction, " is not implemented. Please use harmony, integrated.cca, or integracted.mnn"))
    }
  } else {
    print('Only 1 sample, no need to integrate')
  }
}

# perform dimenstional reduction and clustering
# results are saved 
SeuratDimReduction <- function(se.integrated, dims, group, res = 1, reduction, plots.format){
  se.integrated[["RNA"]] <- JoinLayers(se.integrated[["RNA"]])
  
  
  if ("ATAC" %in% SeuratObject::Assays(se.integrated)) {
    #se.integrated[["ATAC"]] <- JoinLayers(se.integrated[["ATAC"]])
    se.integrated <- FindMultiModalNeighbors(se.integrated, reduction.list = list(reduction, "integrated.lsi"), dims.list = list(1:dims, 2:50))
    #se.integrated <- FindNeighbors(se.integrated, reduction = reduction, dims = 1:dims) # returns KNN graph using the PC or CCA
    se.integrated <- FindClusters(se.integrated, resolution = res, graph.name = "wsnn")
    se.integrated <- RunUMAP(se.integrated, nn.name = "weighted.nn", assay = "RNA")
  } else {
    se.integrated <- FindNeighbors(se.integrated, reduction = reduction, dims = 1:dims) # returns KNN graph using the PC or CCA
    se.integrated <- FindClusters(se.integrated, resolution = res) # find clusters of cells by shared SNN
    se.integrated <- RunUMAP(se.integrated, dims = 1:dims, reduction = reduction) # optimizes the low-dimensional graph representation to be as similar to og graph
  }
  
  p1 <- DimPlot(se.integrated, reduction = "umap", group.by = group, alpha = 0.5) + 
    ggtitle("UMAP with Highlighted Conditions")
  ggsave(paste0("integrated_umap_grouped.", plots.format), p1)
  ggsave(paste0("integrated_umap_grouped.", "svg"), p1)
  
  
  p2 <- DimPlot(se.integrated, reduction = "umap", split.by = group, alpha = 0.5) +
    ggtitle("UMAP Split by Condition and Highlighted by Sample")
  ggsave(paste0("integrated_umap_split.", plots.format), p2)
  ggsave(paste0("integrated_umap_split.", "svg"), p2)

  p3 <- DimPlot(se.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, alpha = 0.5) +
    ggtitle("UMAP Unlabelled Seurat Clusters")
  ggsave(paste0("integrated_umap_unlabelled.", plots.format), p3)
  ggsave(paste0("integrated_umap_unlabelled.", "svg"), p3)

  p4 <- dittoBarPlot(
    object = se.integrated,
    var = "seurat_clusters",
    group.by = "group",
    retain.factor.levels=TRUE,
    main = "Percent of cells for Unlabelled Clusters")
  ggsave(paste0("percent_cells_group_unlabelled.", plots.format), p4)
  ggsave(paste0("percent_cells_group_unlabelled.", "svg"), p4)

  se.integrated <- se.integrated
}

IntegrateAtac <- function(se.merged.preprocessed){
  
  # se.merged.preprocessed[['ATAC']]<- CreateAssay5Object(counts = se.merged.preprocessed@assays[["ATAC"]]@counts,
  #                                                       se.merged.preprocessed@assays[["ATAC"]]@data)

  seurat_obj_list_atac <- SplitObject(se.merged.preprocessed, split.by = "sample")

  seurat_obj_list_atac <- lapply(seurat_obj_list_atac, function(x) {
    DefaultAssay(x) <- "ATAC" # or "RNA" or any other assay name you want to set as default
    return(x)
  })
  
  integration.anchors <- FindIntegrationAnchors(
    object.list = seurat_obj_list_atac,
    reduction = "rlsi",
    dims = 2:50
  )
  print("found anchors")
  # anchor.features = rownames(seurat_obj_list[[1]]), # believe this can be any would this not be needed?

  # integrate LSI embeddings
  DefaultAssay(se.merged.preprocessed) <- "ATAC"
  
  se.merged.preprocessed.atac <- IntegrateEmbeddings(
    anchorset = integration.anchors,
    reductions = se.merged.preprocessed[["lsi"]],
    new.reduction.name = "integrated.lsi",
    dims.to.integrate = 1:50,
    k.weight = FALSE
  )
  print("completed integrated embeddings")
  
  return(se.merged.preprocessed.atac)
}
