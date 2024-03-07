source("./utils/misc.R")

PreprocessingSeurat <- function(seurat_object){
  seurat_object <- NormalizeData(seurat_object) # feature counts for each cell is divided by total counts for the cell; not based on sample
  seurat_object <- FindVariableFeatures(seurat_object) # chooses top variable features (genes); all cells are merged together anyways
  
  # regressing percent.mt for now, for cell cycle could hold importance for cells undergoing differentiation
  # also automating whether the difference in clustering is affected by cell cycle phase requires visual checks
  seurat_object <- ScaleData(seurat_object, vars.to.regress = c("percent.mt")) # scale and center the data, while regressing specified features
  seurat_object <- RunPCA(seurat_object) # dimensional reduction -> large set of genes into smaller metafeature clusters
}

# integratation
IntegrateSamples <- function(seurat_obj_list, group){
  if (length(seurat_obj_list) >= 2){
    se.merged <- merge(seurat_obj_list[[1]], seurat_obj_list[c(2:length(seurat_obj_list))])
    se.merged[["RNA"]] <- split(se.merged[["RNA"]], f = se.merged$group) # whether the data is split by sample or by treatment does not matter for downstream analysis
    se.merged <- UpdateSeuratObject(se.merged)
    
    se.merged.preprocessed <- PreprocessingSeurat(se.merged)
    
    # CCA integration background; https://hbctraining.github.io/scRNA-seq_online/lessons/06_integration.html
    # - It is a form of PCA, in that it identifies the greatest sources of variation in the data, but only if it is shared or conserved across the conditions/groups
    DefaultAssay(se.merged) <- "RNA"
    se.integrated <- IntegrateLayers(object = se.merged.preprocessed, method = CCAIntegration,
                                     orig.reduction = "pca", new.reduction = "integrated.cca",
                                     verbose = T)
  } else {
    print('Only 1 sample, no need to integrate')
  }
}

# perform dimenstional reduction and clustering
# results are saved 
SeuratDimReduction <- function(se.integrated, dims, plot.path, group, res = 1, reduction = "integrated.cca"){
  se.integrated[["RNA"]] <- JoinLayers(se.integrated[["RNA"]])
  se.integrated <- FindNeighbors(se.integrated, reduction = reduction, dims = dims) # returns KNN graph using the PC or CCA
  se.integrated <- FindClusters(se.integrated, resolution = res) # find clusters of cells by shared SNN
  se.integrated <- RunUMAP(se.integrated, dims = dims, reduction = reduction) # optimizes the low-dimensional graph representaiton to be as similar to og graph
  
  DimPlot(se.integrated, reduction = "umap", group.by = group) %>% PrintSave("integrated_umap_grouped.pdf", plot.path)
  DimPlot(se.integrated, reduction = "umap", split.by = group) %>% PrintSave("integrated_umap_split.pdf", plot.path)
  DimPlot(se.integrated, reduction = "umap", group.by = "seurat_clusters") %>% PrintSave("integrated_umap_unlabelled.pdf", plot.path)
  
  se.integrated <- se.integrated
}