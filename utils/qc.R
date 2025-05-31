##### fxns for QC #####

#' Perform ambient RNA correction using soupX over
#' 
#' @param matriceList List that contains the filtered and raw matrices from 
#' Cellranger data
#' @param downsampleNum Denotes the number to downsample the counts matrix
#'
#' @return Folder with the given sample name, with the corrected counts
#' 
AmbientRNARemoval <- function(matriceList, 
                              downsampleNum
                              ) {

  print(paste("Loading ", matriceList[1], " and ", matriceList[2], sep=''))
  
  # check if it is a multiome experiment format
  filt.matrix <- Read10X_h5(matriceList[1],use.names = T)
  raw.matrix <- Read10X_h5(matriceList[2],use.names = T)
  
  if (!is.list(filt.matrix)){
    se <- CreateSeuratObject(counts = filt.matrix) # create seurat object
    
    print("Providing Soup")
    soup.channel <- SoupChannel(raw.matrix, filt.matrix)   # soup channel stpres all info related to a single 10X channel
  } else {
    se <- CreateSeuratObject(counts = filt.matrix$`Gene Expression`) # create seurat object
    
    print("Providing Soup multiome")
    soup.channel <- SoupChannel(raw.matrix$`Gene Expression`, filt.matrix$`Gene Expression`)   # soup channel stpres all info related to a single 10X channel
  }
  
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
  
  if(downsampleNum != 0){
    rand =  sample(1:ncol(adj.matrix), downsampleNum)
    adj.matrix = adj.matrix[, rand]
  }
  
  # should work with cellranger outputs -> just back tracks a set amount of times, so as long as folder struc is correct
  sample.name <- sub(".*\\/(.*)\\/.*\\/.*", "\\1", matriceList[1])
  group.name <- sub(".*\\/(.*)\\/.*\\/.*\\/.*", "\\1", matriceList[1])
  dir.create(paste0("./qc/"))
  dir.create(paste0("./qc/", group.name))
  
  print(paste("Saving under ./qc/", group.name, '/', sample.name, "_soupx",sep=''))
  DropletUtils:::write10xCounts(paste('./qc/', group.name, '/', sample.name, "_soupx",sep=''), adj.matrix, overwrite = TRUE) # name will be the fo;der before the /outs/ folder
}

#' Performs QC filtering
#' 
#' This function will perform automatic filtering using 3 median average 
#' deviations (MADs) for the percentage of mitochnodrial reads, number of UMIs, 
#' and number of genes. Also, if the data is multiomic (scRNA and scATAC), it
#' will filter for the total number of fragments in peaks, TSS enrichment score,
#' and nucleosomal signal.
#' 
#' @param seuratObj A seurat object 
#' @param species Species name, either "musmusculus" or "homosapiens" only for now
#' @param atac Boolean true or false if data is multiomic
#' @param mitoCutoff Set the minimum cutoff used for filtering mitochondrial 
#' percentages
#' @param plotsFormat Extension for the plots
#' 
#' @return Filtered seurat object
#' 
BasicQC <- function(seuratObj, 
                    species, 
                    atac, 
                    mitoCutoff, 
                    plotsFormat
                    ) {
  
  print("Removing low quality cells based on MAD thresholds")
  
  if(all(grepl("^ENS", rownames(seuratObj)))){ #checks if geneIDs are ensembl IDs; need to convert it into gene symbols to filter for mitochondrial genes
    ens.to.symbols <- NULL
    if (species == "musmusculus"){
      library(EnsDb.Mmusculus.v79)
      ens.to.symbols <- as.data.frame(mapIds(EnsDb.Mmusculus.v79, keys = rownames(seuratObj@assays[["RNA"]]@counts),
                                             column = c('SYMBOL'), keytype = 'GENEID'))
    } else if (species == "homosapiens"){
      library(EnsDb.Hsapiens.v86)
      ens.to.symbols <- as.data.frame(mapIds(EnsDb.Hsapiens.v86, keys = rownames(seuratObj@assays[["RNA"]]@counts),
                                             column = c('SYMBOL'), keytype = 'GENEID'))
    }
    names(ens.to.symbols)[1] <- "SYMBOL"
    seuratObj@assays$RNA@meta.features <- merge(seuratObj@assays$RNA@meta.features, ens.to.symbols, by=0, all.x=TRUE) %>% select(1,3) #save under meta.features, allows mapping for later on if needed
    mt.to.calc <- subset(seuratObj@assays$RNA@meta.features, 
                         SYMBOL %in% grep("^mt-", seuratObj@assays$RNA@meta.features$SYMBOL, 
                                          value = TRUE, ignore.case = TRUE))[,1]
    
    seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, features = mt.to.calc)
  } else {
    seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = paste(c("^mt-","^MT-"), collapse="|"))
  }
  
  VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)   # view distribution and to spot any obvious outliers; not saved so can remove
  
  # perform MAD to determine automatic cutoffs
  
  if(atac == "yes"){
    DefaultAssay(seuratObj) <- "ATAC"
    
    seuratObj <- NucleosomeSignal(seuratObj)
    seuratObj <- TSSEnrichment(seuratObj, fast=FALSE)
    
    DefaultAssay(seuratObj) <- "RNA"
  }
  Cell.QC.Stat <- seuratObj@meta.data
  
  ###### Percent mitochondrial filtering #####
  if (mitoCutoff == 0){
    max.mito.thr <- median(Cell.QC.Stat$percent.mt) + 3*mad(Cell.QC.Stat$percent.mt) #looking where to make the cutoff for the max
    if (max.mito.thr < 5){
    	max.mito.thr <- 5
    }
  } else {
    max.mito.thr <- mitoCutoff
  }
  
  p1 <- ggplot(Cell.QC.Stat, aes(x=nFeature_RNA, y=percent.mt)) +
    geom_point(alpha=0.3) +
    geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
    annotate(geom = "text", label = paste0(as.numeric(table(Cell.QC.Stat$percent.mt > max.mito.thr)[2])," cells removed\n",
                                           as.numeric(table(Cell.QC.Stat$percent.mt > max.mito.thr)[1])," cells remain"), x = 6000, y = 0.1) + 
    ggtitle(paste0(seuratObj@misc[[1]], " QC: \nPercentage of Mitochondrial reads per Cell")) +
    xlab("Number of Genes") +
    ylab("Percent of Mitochondrial reads")
  ggsave(paste0(seuratObj@misc[[1]], "_percent_mt.", plotsFormat), plot=p1)
  ggsave(paste0(seuratObj@misc[[1]], "_percent_mt.jpeg"), plot=p1)
  
  Cell.QC.Stat.mt <- filter(Cell.QC.Stat, percent.mt <= max.mito.thr)
  if(nrow(Cell.QC.Stat.mt)/nrow(Cell.QC.Stat) > .5){
    Cell.QC.Stat <- Cell.QC.Stat.mt
  } else {
    stop("Over 50% of your cells were removed because of a mitochondrial cutoff, please check your data")
  }
  
  #ggMarginal(p1, type = "histogram", fill="lightgrey", bins=100) # I think extranaeous but can add later if want
  #min.mito.thr <- median(ifnb$percent.mt) - 3*mad(ifnb$percent.mt) # Many scRNA-seq data doesnt bother with finding a min percent.mt
  
  ##### Filtering by number of genes and number of UMIs #####
  min.Genes.thr <- median(log10(Cell.QC.Stat$nFeature_RNA)) - 3*mad(log10(Cell.QC.Stat$nFeature_RNA))
  max.Genes.thr <- median(log10(Cell.QC.Stat$nFeature_RNA)) + 3*mad(log10(Cell.QC.Stat$nFeature_RNA))
  
  min.nUMI.thr <- median(log10(Cell.QC.Stat$nCount_RNA)) - 3*mad(log10(Cell.QC.Stat$nCount_RNA))
  max.nUMI.thr <- median(log10(Cell.QC.Stat$nCount_RNA)) + 3*mad(log10(Cell.QC.Stat$nCount_RNA))
  
  p2 <- ggplot(Cell.QC.Stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
    geom_point(alpha=0.3) +
    geom_smooth(method="lm") +
    geom_hline(aes(yintercept = min.Genes.thr), colour = "green", linetype = 2) +
    geom_hline(aes(yintercept = max.Genes.thr), colour = "green", linetype = 2) +
    geom_vline(aes(xintercept = min.nUMI.thr), colour = "red", linetype = 2) +
    geom_vline(aes(xintercept = max.nUMI.thr), colour = "red", linetype = 2) +
    annotate(geom = "text", label = paste0(as.numeric(table(log10(Cell.QC.Stat$nFeature_RNA) > max.Genes.thr | 
                                                              log10(Cell.QC.Stat$nFeature_RNA) < min.Genes.thr | 
                                                              log10(Cell.QC.Stat$nCount_RNA) > max.nUMI.thr | 
                                                              log10(Cell.QC.Stat$nCount_RNA) < min.nUMI.thr)[2])," cells removed\n",
                                           as.numeric(table(log10(Cell.QC.Stat$nFeature_RNA) > max.Genes.thr | 
                                                              log10(Cell.QC.Stat$nFeature_RNA) < min.Genes.thr | 
                                                              log10(Cell.QC.Stat$nCount_RNA) > max.nUMI.thr | 
                                                              log10(Cell.QC.Stat$nCount_RNA) < min.nUMI.thr)[1])," cells remain"), x = 3, y = 3) +
    ggtitle(paste0(seuratObj@misc[[1]], " QC: \nNumber of Genes vs Number of UMIs")) +
    xlab("Log10 Transformed Number of UMIs") +
    ylab("Log10 Transformed Number of Genes")
  ggsave(paste0(seuratObj@misc[[1]], "_nGenes_nUMI.", plotsFormat), plot=p2)
  ggsave(paste0(seuratObj@misc[[1]], "_nGenes_nUMI.jpeg"), plot=p2)
  
  if(atac == "yes"){
    DefaultAssay(seuratObj) <- "ATAC"

    # https://stuartlab.org/signac/articles/pbmc_vignette view QC section for QC filtering techniques
    min.peaks.thr <- median(log10(Cell.QC.Stat$nCount_ATAC)) - 3*mad(log10(Cell.QC.Stat$nCount_ATAC))
    max.peaks.thr <- median(log10(Cell.QC.Stat$nCount_ATAC)) + 3*mad(log10(Cell.QC.Stat$nCount_ATAC))
    
    max.nuc.thr <- median(Cell.QC.Stat$nucleosome_signal) + 3*mad(Cell.QC.Stat$nucleosome_signal)
    min.TSS.thr <- median(Cell.QC.Stat$TSS.enrichment) - 3*mad(Cell.QC.Stat$TSS.enrichment)
    
    # TSS score graph
    seuratObj$high.tss <- ifelse(seuratObj$TSS.enrichment > min.TSS.thr, 'High', 'Low')
	unique_tss <- unique(seuratObj$high.tss) 	
   
    if (length(unique_tss) == 1) {
  # If all High or all Low, plot without group.by and add label for that status
  p3 <- TSSPlot(seuratObj) + NoLegend() +
    ggtitle(paste0(seuratObj@misc[[1]], " QC: \nTranscriptional start site (TSS) enrichment score")) +
    labs(tag = paste0("All cells have TSS status: ", unique_tss)) +
    theme(plot.tag.position = c(0, 0))
} else {
  # If both High and Low present, plot with group.by and normal labels
  p3 <- TSSPlot(seuratObj, group.by = 'high.tss') + NoLegend() +
    ggtitle(paste0(seuratObj@misc[[1]], " QC: \nTranscriptional start site (TSS) enrichment score")) +
    labs(tag = paste0(as.numeric(table(Cell.QC.Stat$TSS.enrichment < min.TSS.thr)[2]), " cells removed\n",
                       as.numeric(table(Cell.QC.Stat$TSS.enrichment < min.TSS.thr)[1]), " cells remain")) +
    theme(plot.tag.position = c(0, 0))
}
    ggsave(paste0(seuratObj@misc[[1]], "_tss.", plotsFormat), plot=p3)
    ggsave(paste0(seuratObj@misc[[1]], "_tss.jpeg"), plot=p3)
    
    # Nucleosome signal graph
	# Create the nucleosome_group column
seuratObj$nucleosome_group <- ifelse(seuratObj$nucleosome_signal > max.nuc.thr,
                                     paste0('NS >', max.nuc.thr),
                                     paste0('NS <', max.nuc.thr))

# Identify the unique value if there’s only one
unique_nuc_group <- unique(seuratObj$nucleosome_group)

# Check if there’s only one unique group
if (length(unique_nuc_group) == 1) {
  # If all NS above or below the threshold, plot without group.by and add label
  p4 <- FragmentHistogram(object = seuratObj) +
    ggtitle(paste0(seuratObj@misc[[1]], " QC: \nNucleosome banding pattern")) +
    labs(tag = paste0("All cells have nucleosome group: ", unique_nuc_group)) +
    theme(plot.tag.position = c(0, 0))
} else {
  # If both groups are present, plot with group.by and normal labels
  p4 <- FragmentHistogram(object = seuratObj, group.by = 'nucleosome_group') +
    ggtitle(paste0(seuratObj@misc[[1]], " QC: \nNucleosome banding pattern")) +
    labs(tag = paste0(as.numeric(table(Cell.QC.Stat$nucleosome_signal > max.nuc.thr)[2]), " cells removed\n",
                       as.numeric(table(Cell.QC.Stat$nucleosome_signal > max.nuc.thr)[1]), " cells remain")) +
    theme(plot.tag.position = c(0, 0))
}

    ggsave(paste0(seuratObj@misc[[1]], "_nucleosome_signal.", plotsFormat), plot=p4)
    ggsave(paste0(seuratObj@misc[[1]], "_nucleosome_signal.jpeg"), plot=p4)
    
    seuratObj$log10_nCount_ATAC <- log10(seuratObj$nCount_ATAC)
    p5 <- VlnPlot(object = seuratObj, features = "log10_nCount_ATAC", 
                  pt.size = 1, alpha = 0.2) + 
      geom_hline(yintercept = c(min.peaks.thr, max.peaks.thr), linetype = "dashed", color = "red") + 
      ylab("log10(nCount_ATAC)") +
      ggtitle(paste0(seuratObj@misc[[1]], " QC: \nTotal number of fragments in peaks")) + 
      labs(tag = paste0(as.numeric(table(log10(Cell.QC.Stat$nCount_RNA) > max.peaks.thr | 
                                log10(Cell.QC.Stat$nCount_ATAC) < min.peaks.thr)[2])," cells removed\n",
             as.numeric(table(log10(Cell.QC.Stat$nCount_ATAC) > max.peaks.thr | 
                                log10(Cell.QC.Stat$nCount_ATAC) < min.peaks.thr)[1])," cells remain")) +  
      theme(plot.tag.position = c(0, 0))
    ggsave(paste0(seuratObj@misc[[1]], "_ncount_atac.", plotsFormat), plot=p5) 
    ggsave(paste0(seuratObj@misc[[1]], "_ncount_atac.jpeg"), plot=p5)
    
    Cell.QC.Stat.nfeature.numi <- Cell.QC.Stat %>% # removing the low quality cells from scATAC standards
      filter(log10(nCount_ATAC) >= min.peaks.thr) %>%
      filter(log10(nCount_ATAC) <= max.peaks.thr) %>%
      filter((nucleosome_signal) <= max.nuc.thr) %>%
      filter((TSS.enrichment) >= min.TSS.thr) %>%
      filter(log10(nFeature_RNA) >= min.Genes.thr) %>%
      filter(log10(nFeature_RNA) <= max.Genes.thr) %>%
      filter(log10(nCount_RNA) >= min.nUMI.thr) %>%
      filter(log10(nCount_RNA) <= max.nUMI.thr)
    
    DefaultAssay(seuratObj) <- "RNA"
  } else {
    
    Cell.QC.Stat.nfeature.numi <- Cell.QC.Stat %>% # removing the low quality cells
      filter(log10(nFeature_RNA) >= min.Genes.thr) %>%
      filter(log10(nFeature_RNA) <= max.Genes.thr) %>%
      filter(log10(nCount_RNA) >= min.nUMI.thr) %>%
      filter(log10(nCount_RNA) <= max.nUMI.thr)
  }
  
  if(nrow(Cell.QC.Stat.nfeature.numi)/nrow(Cell.QC.Stat) > .5){
    Cell.QC.Stat <- Cell.QC.Stat.nfeature.numi
  } else {
    stop("Over 50% of your cells were removed because of the number of UMIs and 
         genes cutoff, please check your double check your data")
  }
  
  print("Filtered low-quality cells based of MAD")

  Cell.QC.Stat.bc <- rownames_to_column(Cell.QC.Stat, "barcode")
  filtered.seurat <- subset(seuratObj, cells = Cell.QC.Stat.bc$barcode)   # filtering in Seurat object
}

#' Identify doublets using scDblFinder
#' 
#' @param seuratObj A seurat object  
#' @param atac Boolean true or false if data is multiomic
#'
#' @return Seurat object with metadata denoting singlets or doublets for each cell
#' 
DoubletQC <- function(seuratObj, atac){
  sce <- as.SingleCellExperiment(seuratObj)
  sce <- scDblFinder(sce, clusters=FALSE) # generates random doublets -> generates a new PCA -> creates a kNN network
  # training an iterative classifier on the neighborhood of real cells and artificial doublets

  se <- as.Seurat(sce, counts = "counts", data = NULL)
  
  if (atac == "yes"){
    se[["ATAC"]] <- seuratObj[["ATAC"]]
  }
  return(se)
}
