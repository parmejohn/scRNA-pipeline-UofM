set.seed(333)

#' Perform processing on seurat object before integration
#' 
#' Includes normalizing, finding variable featues, scaling, and running PCA on
#' the data. Please note that SCTransform carries out the same steps, but also 
#' is listed to be more efficient than the method used below.
#' 
#' If the ATAC assay is seen in the object, it will also carry out the
#' processing for steps for the scATAC-seq data.
#' - Term Freq-Inverse Document Frequency (TF-IDF) normalization = normalize cell 
#' sequencing depth, and across peaks to show higher values for rare peaks
#'- Dimensional reduction = Singular Value Decomposition (SVD) on the normalized 
#' matrix for the features selected
#' - similar to PCA output, but for much more sparse data
#'- when performed in this order = latent semantic indexing (LSI)
#' 
#' @param seuratObj A seurat object  
#'
#' @return Seurat object that has undergone basic seurat processing before 
#' integration
#' 
PreprocessingSeurat <- function(seuratObj){
  DefaultAssay(seuratObj) <- "RNA"
  
  seuratObj <- NormalizeData(seuratObj) # feature counts for each cell is divided by total counts for the cell; not based on sample
  seuratObj <- FindVariableFeatures(seuratObj) # chooses top variable features (genes); all cells are merged together anyways
  
  # regressing percent.mt for now, for cell cycle could hold importance for cells undergoing differentiation
  # also automating whether the difference in clustering is affected by cell cycle phase requires visual checks
  seuratObj <- ScaleData(seuratObj, vars.to.regress = c("percent.mt")) # scale and center the data, while regressing specified features
  seuratObj <- RunPCA(seuratObj) # dimensional reduction -> large set of genes into smaller metafeature clusters
  
  if ("ATAC" %in% SeuratObject::Assays(seuratObj)) {
    DefaultAssay(seuratObj) <- "ATAC"
    seuratObj <- RunTFIDF(seuratObj, assay = "ATAC")

    seuratObj <- FindTopFeatures(seuratObj, min.cutoff = "q5", assay = "ATAC")

    seuratObj <- RunSVD(seuratObj, assay = "ATAC")

    DefaultAssay(seuratObj) <- "RNA"
  }
  return(seuratObj)
}

#' Perform integration on a list of seurat objects
#' 
#' Integrates your samples (seurat objects) via harmony, CCA, or FastMNN.
#' If the ATAC assay is found, it will perform multimodal integration to account
#' for the ATAC assays.
#' https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis
#' 
#' @param seuratObjList A seurat object  
#' @param reduction A reduction method to use for integration
#'
#' @return List of seurat objects that have undergone the basic seurat 
#' processing steps (normalized, identified variable features/genes, scaled, and
#' dimensional reduction)
#' 
IntegrateSamples <- function(seuratObjList, 
                             reduction
                             ) {
  if (length(seuratObjList) >= 2){
    se.merged <- merge(seuratObjList[[1]], seuratObjList[c(2:length(seuratObjList))])
    se.merged[["RNA"]] <- split(se.merged[["RNA"]], f = se.merged$group) # whether the data is split by sample or by treatment does not matter for downstream analysis

    se.merged <- UpdateSeuratObject(se.merged)
    se.merged.preprocessed <- PreprocessingSeurat(se.merged)
    
    if ("ATAC" %in% SeuratObject::Assays(se.merged.preprocessed)) {
      se.merged.preprocessed.atac <- IntegrateAtac(se.merged.preprocessed)
      se.merged.preprocessed[['ATAC']] <- se.merged.preprocessed.atac[['ATAC']]
      se.merged.preprocessed[['integrated.lsi']] <- se.merged.preprocessed.atac[['integrated.lsi']]
      # DefaultAssay(se.merged.preprocessed) <- "RNA"
      # se.integrated <- se.merged.preprocessed
    } 
    
    if (reduction == "harmony"){
      DefaultAssay(se.merged.preprocessed) <- "RNA"
      se.integrated <- IntegrateLayers(object = se.merged.preprocessed, method = HarmonyIntegration,
                                       orig.reduction = "pca", new.reduction = "harmony",
                                       assay = "RNA",
                                       verbose = T)
      
    } else if (reduction == "integrated.cca"){
      DefaultAssay(se.merged.preprocessed) <- "RNA"
      se.integrated <- IntegrateLayers(object = se.merged.preprocessed, method = CCAIntegration,
                                       orig.reduction = "pca", new.reduction = "integrated.cca",
                                       assay = "RNA",
                                       verbose = T)
      
    } else if (reduction == "integrated.mnn"){
      DefaultAssay(se.merged.preprocessed) <- "RNA"
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

#' Identifies neighbours, clusters, and performs UMAP dimentional reduction on
#' the reduction method used from integration
#' 
#' @param seuratIntegrated An integrated seurat object
#' @param dims Dimensions to use in the machine learning methods (Clustering, UMAP
#' etc.)
#' @param res Controls the granularity of the clustering. Higher = more clusters
#' @param reduction Name of the reduction used for integration  
#' @param plots.format Extension for the plots
#'
#' @return Creates unlabelled UMAP and stacked bar plots. Also returns a seurat
#' object which has undergone the steps listed above
#' 
SeuratDimReduction <- function(seuratIntegrated, 
                               dims, 
                               group, 
                               res = 1, 
                               reduction, 
                               plots.format
                               ) {
  seuratIntegrated[["RNA"]] <- JoinLayers(seuratIntegrated[["RNA"]])
  
  
  if ("ATAC" %in% SeuratObject::Assays(seuratIntegrated)) {
    #seuratIntegrated[["ATAC"]] <- JoinLayers(seuratIntegrated[["ATAC"]])
    seuratIntegrated <- FindMultiModalNeighbors(seuratIntegrated, reduction.list = list(reduction, "integrated.lsi"), dims.list = list(1:dims, 2:50))
    #seuratIntegrated <- FindNeighbors(seuratIntegrated, reduction = reduction, dims = 1:dims) # returns KNN graph using the PC or CCA
    seuratIntegrated <- FindClusters(seuratIntegrated, resolution = res, graph.name = "wsnn")
    seuratIntegrated <- RunUMAP(seuratIntegrated, nn.name = "weighted.nn", assay = "RNA")
  } else {
    seuratIntegrated <- FindNeighbors(seuratIntegrated, reduction = reduction, dims = 1:dims) # returns KNN graph using the PC or CCA
    seuratIntegrated <- FindClusters(seuratIntegrated, resolution = res) # find clusters of cells by shared SNN
    seuratIntegrated <- RunUMAP(seuratIntegrated, dims = 1:dims, reduction = reduction) # optimizes the low-dimensional graph representation to be as similar to og graph
  }
  
  p1 <- DimPlot(seuratIntegrated, reduction = "umap", group.by = group, alpha = 0.5) + 
    ggtitle("UMAP with Highlighted Conditions")
  ggsave(paste0("integrated_umap_grouped.", plots.format), p1)
  ggsave(paste0("integrated_umap_grouped.", "jpeg"), p1)
  
  
  p2 <- DimPlot(seuratIntegrated, reduction = "umap", split.by = group, alpha = 0.5) +
    ggtitle("UMAP Split by Condition and Highlighted by Sample")
  ggsave(paste0("integrated_umap_split.", plots.format), p2)
  ggsave(paste0("integrated_umap_split.", "jpeg"), p2)

  p3 <- DimPlot(seuratIntegrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, alpha = 0.5) +
    ggtitle("UMAP Unlabelled Seurat Clusters")
  ggsave(paste0("integrated_umap_unlabelled.", plots.format), p3)
  ggsave(paste0("integrated_umap_unlabelled.", "jpeg"), p3)

  p4 <- dittoBarPlot(
    object = seuratIntegrated,
    var = "seurat_clusters",
    group.by = "group",
    retain.factor.levels=TRUE,
    main = "Percent of cells for Unlabelled Clusters")
  ggsave(paste0("percent_cells_group_unlabelled.", plots.format), p4)
  ggsave(paste0("percent_cells_group_unlabelled.", "jpeg"), p4)

  seuratIntegrated <- seuratIntegrated
}

#' Perform LSI integration for ATAC-seq data. This data will be used later on 
#' for multimodaling with the RNA-seq data
#' 
#' @param se.merged.preprocessed A seurat object that has undergone preprocessing
#' specifically for ATAC-seq
#'
#' @return Seurat object with integrated ATAC-seq data
#' 
IntegrateAtac <- function(se.merged.preprocessed){

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

  # integrate LSI embeddings
  DefaultAssay(se.merged.preprocessed) <- "ATAC"
  
  se.merged.preprocessed.atac <- IntegrateEmbeddings(
    anchorset = integration.anchors,
    reductions = se.merged.preprocessed[["lsi"]],
    new.reduction.name = "integrated.lsi",
    dims.to.integrate = 1:50,
    k.weight = FALSE
  )

  return(se.merged.preprocessed.atac)
}
