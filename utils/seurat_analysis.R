#source(paste0(dirname(dirname(dirname(getwd()))),"/utils/misc.R"))
set.seed(333)

PreprocessingSeurat <- function(seurat_object){
  seurat_object <- NormalizeData(seurat_object) # feature counts for each cell is divided by total counts for the cell; not based on sample
  seurat_object <- FindVariableFeatures(seurat_object) # chooses top variable features (genes); all cells are merged together anyways
  
  # regressing percent.mt for now, for cell cycle could hold importance for cells undergoing differentiation
  # also automating whether the difference in clustering is affected by cell cycle phase requires visual checks
  seurat_object <- ScaleData(seurat_object, vars.to.regress = c("percent.mt")) # scale and center the data, while regressing specified features
  seurat_object <- RunPCA(seurat_object) # dimensional reduction -> large set of genes into smaller metafeature clusters
}

# integratation
IntegrateSamples <- function(seurat_obj_list, group, reduction){
  if (length(seurat_obj_list) >= 2){
    se.merged <- merge(seurat_obj_list[[1]], seurat_obj_list[c(2:length(seurat_obj_list))])
    se.merged[["RNA"]] <- split(se.merged[["RNA"]], f = se.merged$group) # whether the data is split by sample or by treatment does not matter for downstream analysis
    se.merged <- UpdateSeuratObject(se.merged)
    
    se.merged.preprocessed <- PreprocessingSeurat(se.merged)
    
    # CCA integration background; https://hbctraining.github.io/scRNA-seq_online/lessons/06_integration.html
    # - It is a form of PCA, in that it identifies the greatest sources of variation in the data, but only if it is shared or conserved across the conditions/groups
    DefaultAssay(se.merged) <- "RNA"
    if (reduction == "harmony"){
      se.integrated <- IntegrateLayers(object = se.merged.preprocessed, method = HarmonyIntegration,
                                       orig.reduction = "pca", new.reduction = "harmony",
                                       verbose = T)
    } else if (reduction == "integrated.cca"){
      se.integrated <- IntegrateLayers(object = se.merged.preprocessed, method = CCAIntegration,
                                       orig.reduction = "pca", new.reduction = "integrated.cca",
                                       verbose = T)
    } else if (reduction == "integrated.mnn"){
      se.integrated <- IntegrateLayers(object = se.merged.preprocessed, method = FastMNNIntegration,
                                       orig.reduction = "pca", new.reduction = "integrated.mnn",
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
SeuratDimReduction <- function(se.integrated, dims, group, res = 1, reduction){
  se.integrated[["RNA"]] <- JoinLayers(se.integrated[["RNA"]])
  se.integrated <- FindNeighbors(se.integrated, reduction = reduction, dims = dims) # returns KNN graph using the PC or CCA
  se.integrated <- FindClusters(se.integrated, resolution = res) # find clusters of cells by shared SNN
  se.integrated <- RunUMAP(se.integrated, dims = dims, reduction = reduction) # optimizes the low-dimensional graph representaiton to be as similar to og graph
  
  p1 <- DimPlot(se.integrated, reduction = "umap", group.by = group, alpha = 0.5) + 
    ggtitle("UMAP with Highlighted Conditions")
  PrintSave(p1, "integrated_umap_grouped.pdf")
  
  p2 <- DimPlot(se.integrated, reduction = "umap", split.by = group, alpha = 0.5) +
    ggtitle("UMAP Split by Condition and Highlighted by Sample")
  PrintSave(p2, "integrated_umap_split.pdf")
  
  p3 <- DimPlot(se.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, alpha = 0.5) +
    ggtitle("UMAP Unlabelled Seurat Clusters")
  PrintSave(p3, "integrated_umap_unlabelled.pdf")
  
  p4 <- dittoBarPlot(
    object = se.integrated,
    var = "seurat_clusters",
    group.by = "group",
    retain.factor.levels=TRUE,
    main = "Percent of cells for Unlabelled Clusters")
  PrintSave(p4, "percent_cells_group_unlabelled.pdf")

  se.integrated <- se.integrated
}

