##### fxns for QC #####
AmbientRNARemoval <- function(pair_list, test){
  #save.loc <- paste(outdir, '/analysis/data/qc/', sep='')
  
  print(paste("Loading ", pair_list[1], "and ", pair_list[2], sep=''))
  filt.matrix <- Read10X_h5(pair_list[1],use.names = T)
  raw.matrix <- Read10X_h5(pair_list[2],use.names = T)
  
  se <- CreateSeuratObject(counts = filt.matrix) # create seurat object
  
  print("Providing Soup")
  soup.channel <- SoupChannel(raw.matrix, filt.matrix)   # soup channel stpres all info related to a single 10X channel
  
  # Quickly perform clustering -> not looking for optimal clustering since this is a filtering step
  # Reports to have better results if some basic clustering is provided, even with basic seurat clustering
  print("Clustering")
  se <- SCTransform(se, verbose = F) # replaces NormalizeData, ScaleData, and FindVariableFeatures
  se <- RunPCA(se, verbose = F) # dimensional reduction
  se <- RunUMAP(se, dims = 1:30, verbose = F)
  se <- FindNeighbors(se, dims = 1:30, verbose = F)
  se <- FindClusters(se, verbose = T)
  
  meta  <- se@meta.data
  umap  <- se@reductions$umap@cell.embeddings
  soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel <- setDR(soup.channel, umap) # used for visualizations if wanted; could remove
  
  print("Writing")
  soup.channel  <- autoEstCont(soup.channel, doPlot=FALSE)   # automatically estimates the contamination fraction
  adj.matrix  <- adjustCounts(soup.channel, roundToInt = T)   # calculate the resulting corrected count matrix with background contamination removed
  
  if(test != 0){
    rand =  sample(1:ncol(adj.matrix), test)
    adj.matrix = adj.matrix[, rand]
  }
  
  # should work with cellranger outputs -> just back tracks a set amount of times, so as long as folder struc is correct
  sample.name <- sub(".*\\/(.*)\\/.*\\/.*", "\\1", pair_list[1])
  group.name <- sub(".*\\/(.*)\\/.*\\/.*\\/.*", "\\1", pair_list[1])
  dir.create(paste0("./qc/"))
  dir.create(paste0("./qc/", group.name))
  
  print(paste("Saving under ./qc/", group.name, '/', sample.name, "_soupx",sep=''))
  DropletUtils:::write10xCounts(paste('./qc/', group.name, '/', sample.name, "_soupx",sep=''), adj.matrix, overwrite = TRUE) # name will be the fo;der before the /outs/ folder
  #DropletUtils:::write10xCounts(paste(save.loc, group.name, '/', sample.name, "_soupx",sep=''), adj.matrix, overwrite = TRUE) # name will be the fo;der before the /outs/ folder
  #print(paste("Saved under ", save.loc, group.name, '/', sample.name, "_soupx",sep=''))
}

BasicQC <- function(seurat_obj){
  
  print("Removing low quality cells based on MAD thresholds")
  
  
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = paste(c("^mt-","^MT-"), collapse="|"))
  
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)   # view distribution and to spot any obvious outliers; not saved so can remove
  
  # perform MAD to determine automatic cutoffs
  Cell.QC.Stat <- seurat_obj@meta.data
  
  ###### Percent mitochondrial filtering #####
  max.mito.thr <- median(Cell.QC.Stat$percent.mt) + 3*mad(Cell.QC.Stat$percent.mt) #looking where to make the cutoff for the max
  
  p1 <- ggplot(Cell.QC.Stat, aes(x=nFeature_RNA, y=percent.mt)) +
    geom_point() +
    geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
    annotate(geom = "text", label = paste0(as.numeric(table(Cell.QC.Stat$percent.mt > max.mito.thr)[2])," cells removed\n",
                                           as.numeric(table(Cell.QC.Stat$percent.mt > max.mito.thr)[1])," cells remain"), x = 6000, y = 0.1)
  
  # dir.create(paste0(plot.path, "/qc/"), showWarnings = FALSE)
  # dir.create(paste0(plot.path, "/qc/", seurat_obj@misc[[1]]), showWarnings = FALSE)
  ggsave(paste0(seurat_obj@misc[[1]], "_percent_mt.pdf"), plot=p1)
  
  Cell.QC.Stat <- Cell.QC.Stat %>% filter(percent.mt < max.mito.thr)
  
  #ggMarginal(p1, type = "histogram", fill="lightgrey", bins=100) # I think extranaeous but can add later if want
  #min.mito.thr <- median(ifnb$percent.mt) - 3*mad(ifnb$percent.mt) # Many scRNA-seq data doesnt bother with finding a min percent.mt
  
  ##### Filtering by number of genes and number of UMIs #####
  min.Genes.thr <- median(log10(Cell.QC.Stat$nFeature_RNA)) - 3*mad(log10(Cell.QC.Stat$nFeature_RNA))
  max.Genes.thr <- median(log10(Cell.QC.Stat$nFeature_RNA)) + 3*mad(log10(Cell.QC.Stat$nFeature_RNA))
  
  min.nUMI.thr <- median(log10(Cell.QC.Stat$nCount_RNA)) - 3*mad(log10(Cell.QC.Stat$nCount_RNA))
  max.nUMI.thr <- median(log10(Cell.QC.Stat$nCount_RNA)) + 3*mad(log10(Cell.QC.Stat$nCount_RNA))
  
  p2 <- ggplot(Cell.QC.Stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
    geom_point() +
    geom_smooth(method="lm") +
    geom_hline(aes(yintercept = min.Genes.thr), colour = "green", linetype = 2) +
    geom_hline(aes(yintercept = max.Genes.thr), colour = "green", linetype = 2) +
    geom_vline(aes(xintercept = min.nUMI.thr), colour = "red", linetype = 2) +
    geom_vline(aes(xintercept = max.nUMI.thr), colour = "red", linetype = 2)
  ggsave(paste0(seurat_obj@misc[[1]], "_nGenes_nUMI.pdf"), plot=p2)
  
  
  Cell.QC.Stat <- Cell.QC.Stat %>% # removing the low quality cells
    filter(log10(nFeature_RNA) > min.Genes.thr) %>%
    filter(log10(nFeature_RNA) < max.Genes.thr) %>%
    filter(log10(nCount_RNA) > min.nUMI.thr) %>%
    filter(log10(nCount_RNA) < max.nUMI.thr)
  
  print("Filtered low-quality cells based of MAD")
  #print(paste("QC plots saved under ", seurat_obj@misc[[1]],sep=''))
  
  Cell.QC.Stat.bc <- rownames_to_column(Cell.QC.Stat, "barcode")
  filtered.seurat <- subset(seurat_obj, cells = Cell.QC.Stat.bc$barcode)   # filtering in Seurat object
}

DoubletQC <- function(seurat_obj){
  
  sce <- as.SingleCellExperiment(seurat_obj)
  sce <- scDblFinder(sce, clusters=FALSE) # generates random doublets -> generates a new PCA -> creates a kNN network
  # training an iterative classifier on the neighborhood of real cells and artificial doublets
  se <- as.Seurat(sce, counts = "counts", data = NULL)
  se.singlet <- subset(se, subset = scDblFinder.class  == "singlet") # remove doublets from seurat object
  
}
