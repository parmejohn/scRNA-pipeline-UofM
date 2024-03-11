source("./utils/misc.R")

CellProportionByCluster <- function(se.integrated, plot.path){
  conditions <- unique(se.integrated@meta.data[["group"]])
  group.num <- data.frame(.=character(),
                          Freq=character(), 
                          group=character(), 
                          stringsAsFactors=FALSE)
  
  for(i in conditions){
    cond.num <- subset(se.integrated, group == i) %>% Idents() %>% table() %>% 
      as.data.frame() 
    cond.num$group <- i
    group.num <- rbind(group.num, cond.num)
    # treat.num <- subset(se.integrated, group == ident.2) %>% Idents() %>% table() %>% 
    #   as.data.frame()
    # treat.num$group <- ident.2
  }
  
  # group.num <- rbind(ctrl.num, treat.num) %>% 
  group.num <- dplyr::rename(group.num, cluster = '.')
  
  cond.prop.plot <- group.num %>%
    group_by(cluster) %>%
    mutate(prop = Freq / sum(Freq)) %>%
    ggplot(aes(x = cluster, y = prop)) +
    geom_col(aes(color = group, fill = group), position = position_dodge(0.8), width = 0.7) + 
    theme(axis.text.x = element_text(angle = 90)) +
    ylim(0, 1)
  PrintSave(cond.prop.plot, 'condition_proportion_per_cluster.pdf', plot.path)
}

DifferentialAbundanceMilo <- function(se.integrated, sample, condition, k, d, plot.path, reduced.dims, prop = 0.15){

  DefaultAssay(se.integrated) <- "RNA"
  se.integrated$da.clusters <- Idents(se.integrated)
  
  # make sure that reduced dims is carried over to save on time
  sc.integrated <- as.SingleCellExperiment(se.integrated,  assay = "RNA") # need to convert Seurat to SCE object
  sc.integrated.milo <- Milo(sc.integrated) # make a Milo obj

  # Constructing graph -> computes a k-nearest neighbour graph; Graph vertex = single cell, Edges = neighbors
  # built from PCA -> so if using integrated dataset the PCA will be corrected for using MNN
  sc.integrated.milo.traj <- buildGraph(sc.integrated.milo, k = k, d = d, reduced.dim = reduced.dims) # d = use the corrected reduction; k is increased or decreased according to plotNhoodSizeHist

  # Defining neigborhoods -> group of cells connected by an edge in KNN graph to an index cell
  sc.integrated.milo.traj <- makeNhoods(sc.integrated.milo.traj, prop = prop, k = k, d=d, reduced_dims = reduced.dims, refined = TRUE)
  # If distribution peak is not ideal, have to change previous k values
  # can calc bins automatically
  # how many cells form each neighborhood
  plotNhoodSizeHist(sc.integrated.milo.traj) # official vignette says it wants the distribution to peak between 50-100 but others say 5 x N samples. Using official vignette's suggestion for now

  # counting the cells -> variation in cell numbers between replicates for the same cond to test for DA
  sc.integrated.milo.traj <- countCells(sc.integrated.milo.traj, meta.data = as.data.frame(colData(sc.integrated.milo.traj)), sample=sample)

  # Setting up the design
  sc.integrated.milo.traj.design <- data.frame(colData(sc.integrated.milo.traj))[,c(sample, condition)]
  sc.integrated.milo.traj.design <- distinct(sc.integrated.milo.traj.design)
  rownames(sc.integrated.milo.traj.design) <- sc.integrated.milo.traj.design$sample
  ## Reorder rownames to match columns of nhoodCounts(milo)
  sc.integrated.milo.traj.design <- sc.integrated.milo.traj.design[colnames(nhoodCounts(sc.integrated.milo.traj)), , drop=FALSE]

  # store distances -> store distances btw nearest neighbors -> Milo uses this downstream for FDR correction
  # this step takes long to compute
  print("starting long compute time")
  sc.integrated.milo.traj <- calcNhoodDistance(sc.integrated.milo.traj, d=d, reduced.dim = reduced.dims)
  print("end of long compute time")
  
  da.results <- testNhoods(sc.integrated.milo.traj, design = ~ group, design.df = sc.integrated.milo.traj.design, reduced.dim=reduced.dims)
  
  p1 <- ggplot(da.results, aes(PValue)) + geom_histogram(bins=50) # should be anti-conservative distribution
  p2 <- ggplot(da.results, aes(logFC, -log10(SpatialFDR))) + #volcano plot to see how many neighborhoods are above a fdr threshold
    geom_point() +
    geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)
  
  sc.integrated.milo.traj <- buildNhoodGraph(sc.integrated.milo.traj)
  
  # finding DEGs in DA neighborhoods
  da.results <- annotateNhoods(sc.integrated.milo.traj, da.results, coldata_col = "da.clusters")
  logcounts(sc.integrated.milo.traj) <- log1p(counts(sc.integrated.milo.traj))

  print("Finding DEGs for DA neighborhoods, this may take a while")
  for (i in levels(droplevels(se.integrated@meta.data[["da.clusters"]]))){
    print(paste0("Working on DA DE heatmap for cluster ", i))
    dge_smp <- findNhoodMarkers(sc.integrated.milo.traj, da.results,
                                assay = "counts", gene.offset = FALSE, da.fdr = 0.1,
                                aggregate.samples = TRUE, sample_col = "sample",
                                subset.nhoods = eval(parse(text=paste0("da.results", "$", "da.clusters"))) %in% c(i) # seeing if this paste method works
    )
    dge.smp.filt <- dge_smp %>%
      filter_at(vars(starts_with("adj.P.Val_")), any_vars(. <= 0.01))
    
    markers <- dge.smp.filt[, "GeneID"]
    sc.integrated.milo.traj <- calcNhoodExpression(sc.integrated.milo.traj, subset.row=markers)
    
    da.results.filt <- filter(da.results, da.clusters == i & SpatialFDR < 0.1) # Error in hclust(dist(expr_mat)) : must have n >= 2 objects to cluster -> filtered sets mustve had nothing
    
    if (dim(da.results.filt)[1] >= 2 & length(markers) >= 2){ # have to check if there are any that meet the spatialFDR or marker cutoff to begin with
      p5 <- plotNhoodExpressionDA(sc.integrated.milo.traj, da.results.filt, features = markers,
                            subset.nhoods = eval(parse(text=paste0("da.results.filt", "$", "da.clusters"))) %in% c(i),
                            assay="logcounts",
                            scale_to_1 = TRUE, cluster_features = TRUE, show_rownames = FALSE
      )
      PrintSave(p5, paste0("milo_DA_DE_heatmap_", i, ".pdf"), paste(plot.path, '/da/', sep=''))
    }
  }
  # plots the Differential abundance -> uses UMAP reduction from seurat object
  # - nodes = neighborhoods
  #   - color = fold change
  #   - node layout -> based on index node in UMAP of single cells
  # - edges = number cells shared btw adjacent neighborhoods
  p3 <- plotNhoodGraphDA(sc.integrated.milo.traj, da.results, layout="UMAP",alpha=0.1) # for my UMAP negative FC denoted CTRL and positive values denoted TREAT -> will have to look at your original UMAP for density of cells
  # also can look at the your design -> should be ordered correspondingly; https://github.com/MarioniLab/miloR/issues/81
  
  
  da.results <- annotateNhoods(sc.integrated.milo.traj, da.results, coldata_col = "da.clusters")
  p4 <- plotDAbeeswarm(da.results, group.by = "da.clusters") # gives a distribution view instead of UMAP
  
  dir.create(paste(plot.path, "da", sep=""))
  PrintSave(p1, "milo_pval_distribution.pdf", paste(plot.path, '/da/', sep=''))
  PrintSave(p2, "milo_volcano_plot.pdf", paste(plot.path, '/da/', sep=''))
  PrintSave(p3, "milo_DA_umap.pdf", paste(plot.path, '/da/', sep=''))
  PrintSave(p4, "milo_DA_fc_distribution.pdf", paste(plot.path, '/da/', sep=''))
  #sc.integrated.milo.traj
  se.integrated$da.clusters <- NULL
}
