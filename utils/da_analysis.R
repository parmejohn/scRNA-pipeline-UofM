set.seed(333)

CellProportionByCluster <- function(se.integrated) {
  conditions <- unique(se.integrated@meta.data[["group"]])
  group.num <- data.frame(
    . = character(),
    Freq = character(),
    group = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in conditions) {
    cond.num <-
      subset(se.integrated, group == i) %>% Idents() %>% table() %>%
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
    geom_col(aes(color = group, fill = group),
             position = position_dodge(0.8),
             width = 0.7) +
    theme(axis.text.x = element_text(angle = 90)) +
    ylim(0, 1)
  PrintSave(cond.prop.plot, 'condition_proportion_per_cluster.pdf')
}

#se.integrated <- readRDS("/home/projects/sc_pipelines/scrna_deanne_harmony_low_res/pipeline/analysis/data/se_integrated_auto_label.rds")
# k = 16
# d = 50
# sample = "sample"
# condition = "group"
# reduced.dims = "INTEGRATED.CCA"
# prop = 0.05
# 
# sc.integrated.milo.traj <- readRDS("/home/projects/sc_pipelines/scrna_deanne_harmony_low_res/pipeline/analysis/data/sc_integrated_milo_traj.rds")

DifferentialAbundanceMilo <-
  function(se.integrated,
           sample,
           k,
           d,
           reduced.dims,
           prop = 0.05) {
    DefaultAssay(se.integrated) <- "RNA"
    se.integrated$da.clusters <- Idents(se.integrated)
    #se.integrated@reductions[["harmony"]]@stdev <- as.numeric(apply(se.integrated@reductions[["harmony"]]@cell.embeddings, 2, stats::sd))
    
    # make sure that reduced dims is carried over to save on time
    sc.integrated <-
      as.SingleCellExperiment(se.integrated,  assay = "RNA") # need to convert Seurat to SCE object
    sc.integrated.milo <- Milo(sc.integrated) # make a Milo obj
    
    # Constructing graph -> computes a k-nearest neighbour graph; Graph vertex = single cell, Edges = neighbors
    # built from PCA -> so if using integrated dataset the PCA will be corrected for using MNN
    # sc.integrated.milo.traj <-
    #   buildGraph(
    #     sc.integrated.milo,
    #     k = k,
    #     d = d,
    #     reduced.dim = reduced.dims
    #   ) # d = use the corrected reduction; k is increased or decreased according to plotNhoodSizeHist but usually the same as building knn
    # 
    # # Defining neigborhoods -> group of cells connected by an edge in KNN graph to an index cell
    # sc.integrated.milo.traj <-
    #   makeNhoods(
    #     sc.integrated.milo.traj,
    #     prop = prop,
    #     k = k,
    #     d = d,
    #     reduced_dims = reduced.dims,
    #     refined = TRUE
    #   )
    # # If distribution peak is not ideal, have to change previous k values
    # # can calc bins automatically
    # # how many cells form each neighborhood
    # plotNhoodSizeHist(sc.integrated.milo.traj) # official vignette says it wants the distribution to peak between 50-100 but others say 5 x N samples. Using official vignette's suggestion for now
    # 
    # # counting the cells -> variation in cell numbers between replicates for the same cond to test for DA
    # sc.integrated.milo.traj <-
    #   countCells(sc.integrated.milo.traj,
    #              meta.data = as.data.frame(colData(sc.integrated.milo.traj)),
    #              sample = sample)
    # 
    # # Setting up the design
    # sc.integrated.milo.traj.design <-
    #   data.frame(colData(sc.integrated.milo.traj))[, c(sample, condition)]
    # sc.integrated.milo.traj.design <-
    #   distinct(sc.integrated.milo.traj.design)
    # rownames(sc.integrated.milo.traj.design) <-
    #   sc.integrated.milo.traj.design$sample
    # ## Reorder rownames to match columns of nhoodCounts(milo)
    # sc.integrated.milo.traj.design <-
    #   sc.integrated.milo.traj.design[colnames(nhoodCounts(sc.integrated.milo.traj)), , drop =
    #                                    FALSE]
    # 
    # # store distances -> store distances btw nearest neighbors -> Milo uses this downstream for FDR correction
    # # this step takes long to compute
    # print("starting long compute time")
    # sc.integrated.milo.traj <-
    #   calcNhoodDistance(sc.integrated.milo.traj,
    #                     d = d,
    #                     reduced.dim = reduced.dims)
    # print("end of long compute time")
    # 
    # da.results <-
    #   testNhoods(
    #     sc.integrated.milo.traj,
    #     design = ~ group,
    #     design.df = sc.integrated.milo.traj.design,
    #     reduced.dim = reduced.dims
    #   )
    # 
    ##### testing different method
    ## Build KNN graph neighbourhoods
    milo.obj <- buildGraph(sc.integrated.milo, k=k, d=d, reduced.dim = reduced.dims)
    milo.obj <- makeNhoods(milo.obj, k=k, d=d, refined=T, prop=0.2, refinement_scheme="graph", reduced_dims = reduced.dims)

    ## Count cells in nhoods
    milo.obj <- countCells(milo.obj, samples=sample, meta.data=as.data.frame(colData(milo.obj)))

    ## get the list of conditions
    condition.list <- list()
    condition.list <- append(condition.list, "group")
    if (length(se.integrated@misc[[1]]) > 0 ){
      condition.list <- append(condition.list, se.integrated@misc$co.conditions)
    }
    
    ## Perform miloR analysis on each condition; still integrated and normalized by the main condition however
    for (condition in condition.list){
      
      #colData(milo.obj)[condition] <- relevel(mixedsort(colData(milo.obj)[condition], decreasing=TRUE))
      milo.obj@colData@listData[[condition]] <- factor(milo.obj@colData@listData[[condition]], 
                                                       levels = mixedsort(unique(milo.obj@colData@listData[[condition]]), decreasing=FALSE))
      
      ## Setting up the design
      sc.integrated.milo.traj.design <-
        data.frame(colData(milo.obj))[, c(sample, condition)]
      sc.integrated.milo.traj.design <-
        distinct(sc.integrated.milo.traj.design)
      rownames(sc.integrated.milo.traj.design) <-
        sc.integrated.milo.traj.design$sample
      ## Reorder rownames to match columns of nhoodCounts(milo)
      sc.integrated.milo.traj.design <-
        sc.integrated.milo.traj.design[colnames(nhoodCounts(milo.obj)), , drop =
                                         FALSE]
      
      ## Test for differential abundance
      milo_res <- testNhoods(milo.obj, design=~eval(parse(text=condition)), design.df=sc.integrated.milo.traj.design, fdr.weighting="graph-overlap", reduced.dim = reduced.dims)
      da.results <- milo_res
      write.table(da.results, paste0("da_diff_test_", condition, ".txt"), quote = FALSE,row.names = T, sep = "\t", col.names = T)
      
      sc.integrated.milo.traj <- milo.obj
      
      p1 <-
        ggplot(da.results, aes(PValue)) + 
        geom_histogram(bins = 50) + # should be anti-conservative distribution
        ggtitle(paste0("DA Test: Distribution of p-values for ", condition))
      p2 <-
        ggplot(da.results, aes(logFC,-log10(SpatialFDR))) + #volcano plot to see how many neighborhoods are above a fdr threshold
        geom_point() +
        geom_hline(yintercept = 1) +## Mark significance threshold (10% FDR)
        ggtitle(paste0("Significant Neighborhoods over ", condition))
      PrintSave(p1, paste0("milo_pval_distribution_", condition, ".pdf"))
      PrintSave(p2, paste0("milo_volcano_plot_", condition, ".pdf"))
      
      sc.integrated.milo.traj <-
        buildNhoodGraph(sc.integrated.milo.traj)
      
      # finding DEGs in DA neighborhoods
      da.results <-
        annotateNhoods(sc.integrated.milo.traj, da.results, coldata_col = "da.clusters")
      write.table(da.results, paste0("da_diff_test_", condition,".txt"), quote = FALSE,row.names = T, sep = "\t", col.names = T)
      logcounts(sc.integrated.milo.traj) <-
        log1p(counts(sc.integrated.milo.traj))
      print(paste0("buzzbuzz ", condition))
      #da.results$da.clusters <- as.numeric(da.results$da.clusters) # if default seurat_clusters?
      
      p4 <- plotDAbeeswarm_fixed(da.results, group.by = "da.clusters") # gives a distribution view instead of UMAP
      if (is.ggplot(p4)) {
      	p4 <- p4 + ggtitle(paste0("DA FC Distribution: ", condition))
            PrintSave(p4, paste0("milo_DA_fc_distribution_", condition, ".pdf"))
      }
      
      # print("Finding DEGs for DA neighborhoods, this may take a while") https://marionilab.github.io/miloR/articles/milo_gastrulation.html
      for (i in levels(droplevels(se.integrated@meta.data[["da.clusters"]]))) {
        print(paste0("Working on DA DE heatmap for cluster ", i, "for condition ", condition))
        logcounts(sc.integrated.milo.traj) <-
          log1p(counts(sc.integrated.milo.traj))
        
        n.da <- sum(da.results$SpatialFDR < 0.1)
        dge_smp <- 0
        if (!is.na(n.da) & n.da == 0) {
          dge_smp <- NULL
        } else {
          # dge_smp <- findNhoodMarkers( ### METHOD IS BEING DEPRECATED
          #   sc.integrated.milo.traj,
          #   da.results,
          #   assay = "logcounts",
          #   gene.offset = FALSE,
          #   da.fdr = 0.1,
          #   aggregate.samples = TRUE,
          #   sample_col = sample,
          #   subset.nhoods = da.results$da.clusters %in% c(i)
          # ) # seeing if this paste method works
          dge_smp <- groupNhoods(
            sc.integrated.milo.traj,
            da.results,
            da.fdr = 0.1,
            subset.nhoods = da.results$da.clusters %in% c(i)
          )
          dge_smp <- findNhoodGroupMarkers(
            sc.integrated.milo.traj,
            dge_smp,
            assay = "logcounts",
            gene.offset = FALSE,
            aggregate.samples = TRUE,
            sample_col = sample,
            subset.nhoods = da.results$da.clusters %in% c(i)
          )
        }
        
        if (!is.null(dge_smp)) {
          #paste("findNhoodMarkers")
          dge.smp.filt <- dge_smp %>%
            filter_at(vars(starts_with("adj.P.Val_")), any_vars(. <= 0.01))
          #paste("filt adjpval")
          print("renaming genes")
          markers <- dge.smp.filt[, "GeneID"]
          #markers <- sub("HLA.", "HLA-", markers) #temporary solution
          print("renamed")
          
          logcounts(sc.integrated.milo.traj) <-
            log1p(counts(sc.integrated.milo.traj))
          print("get log counts again")
          sc.integrated.milo.traj <-
            calcNhoodExpression(sc.integrated.milo.traj, subset.row = markers)
          #paste("traj")
          
          da.results.filt <-
            filter(da.results, da.clusters == i & SpatialFDR < 0.1) # Error in hclust(dist(expr_mat)) : must have n >= 2 objects to cluster -> filtered sets mustve had nothing
          #print("dafilt")
          
          if (!is.null(dim(da.results.filt)[1]) & !is.null(markers)) {
            if (dim(da.results.filt)[1] >= 2 &
                length(markers) >= 2) {
              # have to check if there are any that meet the spatialFDR or marker cutoff to begin with
              print("plotting deg")
              p5 <-
                plotNhoodExpressionDA_fixed(
                  sc.integrated.milo.traj,
                  da.results.filt,
                  features = markers,
                  subset.nhoods = da.results$da.clusters %in% c(i),
                  assay = "logcounts",
                  scale_to_1 = TRUE,
                  cluster_features = TRUE,
                  show_rownames = FALSE,
                  group = c(i),
                  condition = condition
                ) + 
                ggtitle(paste0(i, ": DA DEGs for condition ", condition))
              PrintSave(p5, paste0("milo_DA_DE_heatmap_", i, "_", condition, ".pdf"))
            }
          }
        }
      }
      
      # plots the Differential abundance -> uses UMAP reduction from seurat object
      # - nodes = neighborhoods
      #   - color = fold change
      #   - node layout -> based on index node in UMAP of single cells
      # - edges = number cells shared btw adjacent neighborhoods
      p3 <-
        plotNhoodGraphDA(sc.integrated.milo.traj,
                         da.results,
                         layout = "UMAP",
                         alpha = 0.1) +  # for my UMAP negative FC denoted CTRL and positive values denoted TREAT -> will have to look at your original UMAP for density of cells
      # also can look at the your design -> should be ordered correspondingly; https://github.com/MarioniLab/miloR/issues/81
        ggtitle(paste0("DA Analysis UMAP for ", condition))
      
      da.results <-
        annotateNhoods(sc.integrated.milo.traj, da.results, coldata_col = "da.clusters")
      #p4 <- plotDAbeeswarm(da.results, group.by = "da.clusters") # gives a distribution view instead of UMAP
      
      #dir.create(paste(plot.path, "da", sep=""))
      #PrintSave(p1, "milo_pval_distribution.pdf")
      #PrintSave(p2, "milo_volcano_plot.pdf")
      PrintSave(p3, paste0("milo_DA_umap_", condition,".pdf"))
      #PrintSave(p4, "milo_DA_fc_distribution.pdf")
      #sc.integrated.milo.traj
      saveRDS(sc.integrated.milo.traj, paste0("sc_integrated_milo_traj_", condition, ".rds"))
    }
    se.integrated$da.clusters <- NULL
  }
