set.seed(333)

DifferentialAbundanceMilo <-
  function(se.integrated,
           sample,
           k,
           d,
           reduced.dims,
           prop = 0.05,
           fdr.cutoff = 0.05,
           species,
           plots.format) {
    DefaultAssay(se.integrated) <- "RNA"
    se.integrated$da.clusters <- Idents(se.integrated)

    # make sure that reduced dims is carried over to save on time
    sc.integrated <-
      as.SingleCellExperiment(se.integrated,  assay = "RNA") # need to convert Seurat to SCE object
    sc.integrated.milo <- Milo(sc.integrated) # make a Milo obj
    m_df <- msigdbr(species = species, category = "C5", subcategory = "BP") # dont need to reload the dataset every time for GSEA
    fgsea.sets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
    
    # Constructing graph -> computes a k-nearest neighbour graph; Graph vertex = single cell, Edges = neighbors
    # built from PCA -> so if using integrated dataset the PCA will be corrected for using MNN

    # # Defining neigborhoods -> group of cells connected by an edge in KNN graph to an index cell

    # # If distribution peak is not ideal, have to change previous k values
    # # can calc bins automatically
    # # how many cells form each neighborhood
    # plotNhoodSizeHist(sc.integrated.milo.traj) # official vignette says it wants the distribution to peak between 50-100 but others say 5 x N samples. Using official vignette's suggestion for now
    # 
    # # counting the cells -> variation in cell numbers between replicates for the same cond to test for DA
    
    # # store distances -> store distances btw nearest neighbors -> Milo uses this downstream for FDR correction
 
    ##### testing different method
    ## Build KNN graph neighbourhoods
    milo.obj <- buildGraph(sc.integrated.milo, k=k, d=d, reduced.dim = reduced.dims)
    milo.obj <- makeNhoods(milo.obj, k=k, d=d, refined=T, prop=prop, refinement_scheme="graph", reduced_dims = reduced.dims)

    ## Count cells in nhoods
    milo.obj <- countCells(milo.obj, samples=sample, meta.data=as.data.frame(colData(milo.obj)))

    ## get the list of conditions
    condition.list <- list()
    condition.list <- append(condition.list, "group")
    if (length(se.integrated@misc) != 0 ){
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
        geom_hline(yintercept = -log10(fdr.cutoff)) + ## Mark significance threshold (10% FDR)
        ggtitle(paste0("Significant Neighborhoods over ", condition))
      PrintSaveAndSVG(p1, paste0("milo_pval_distribution_", condition), plots.format)
      PrintSaveAndSVG(p2, paste0("milo_volcano_plot_", condition), plots.format)
      
      sc.integrated.milo.traj <-
        buildNhoodGraph(sc.integrated.milo.traj)
      
      # finding DEGs in DA neighborhoods
      da.results <-
        annotateNhoods(sc.integrated.milo.traj, da.results, coldata_col = "da.clusters")
      write.table(da.results, paste0("da_diff_test_", condition,".txt"), quote = FALSE,row.names = T, sep = "\t", col.names = T)
      logcounts(sc.integrated.milo.traj) <-
        log1p(counts(sc.integrated.milo.traj))
      print(paste0("buzzbuzz ", condition))

      p4 <- plotDAbeeswarm_fixed(da.results, group.by = "da.clusters", alpha = fdr.cutoff) # gives a distribution view instead of UMAP
      if (is.ggplot(p4)) {
        p4 <- p4 + ggtitle(paste0("DA FC Distribution: ", condition))
        PrintSaveAndSVG(p4, paste0("milo_DA_fc_distribution_", condition), plots.format)
      }
      
      # print("Finding DEGs for DA neighborhoods, this may take a while") https://marionilab.github.io/miloR/articles/milo_gastrulation.html
      for (i in levels(droplevels(se.integrated@meta.data[["da.clusters"]]))) {
        print(paste0("Working on DA DE heatmap for cluster ", i, " for condition ", condition))
        logcounts(sc.integrated.milo.traj) <-
          log1p(counts(sc.integrated.milo.traj))
        
        n.da <- sum(da.results$SpatialFDR < fdr.cutoff)
        dge_smp <- 0
        if (!is.na(n.da) & n.da == 0) {
          dge_smp <- NULL
        } else if (!(i %in% da.results$da.clusters)){
          dge_smp <- NULL
        } else if (sum(i == da.results$da.clusters) < 2){
          dge_smp <- NULL
        } else {
          dge_smp <- groupNhoods(
            sc.integrated.milo.traj,
            da.results,
            da.fdr = fdr.cutoff,
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
          print("found group markers")
        }
        
        if (!is.null(dge_smp) | length(dge_smp) != 0) {
          dge.smp.filt <- dynamic_filter_function(dge_smp)
          
          ### perform GSEA analyses given the marker genes
          dge.smp.filt.avg.fc <- column_to_rownames(dge.smp.filt, "GeneID")
          dge.smp.filt.avg.fc <- select(dge.smp.filt.avg.fc, c(contains("logFC_")))
          dge.smp.filt.avg.fc$avg_logFC <- rowMeans(dge.smp.filt.avg.fc)
          dge.smp.filt.avg.fc <- rownames_to_column(dge.smp.filt.avg.fc, "gene")
          write.table(dge.smp.filt.avg.fc, paste0("da_", i, "_markers_avg_logfc_", condition,".txt"), quote = FALSE,row.names = T, sep = "\t", col.names = T)
          
          DAGseaComparison(dge.smp.filt.avg.fc, i, condition, fgsea.sets, plots.format)
          
          print("renaming genes")
          markers <- dge.smp.filt[, "GeneID"]
          if (!is.null(markers)){
            print("renamed")
            
            logcounts(sc.integrated.milo.traj) <-
              log1p(counts(sc.integrated.milo.traj))
            print("get log counts again")
            sc.integrated.milo.traj <-
              calcNhoodExpression(sc.integrated.milo.traj, subset.row = markers)
            #paste("traj")
            
            da.results.filt <-
              filter(da.results, da.clusters == i & SpatialFDR < fdr.cutoff) # Error in hclust(dist(expr_mat)) : must have n >= 2 objects to cluster -> filtered sets mustve had nothing
            #print("dafilt")
            
            if (!is.null(dim(da.results.filt)[1]) & !is.null(markers)) {
              if (dim(da.results.filt)[1] >= 2 &
                  length(markers) >= 2) {
                # have to check if there are any that meet the spatialFDR or marker cutoff to begin with
                print("plotting deg")
                #pdf(paste0("milo_DA_DE_heatmap_", i, "_", condition, ".pdf"), width = 8, height = 6)
                  plotNhoodExpressionDA_fixed(
                    x = sc.integrated.milo.traj,
                    da.res = da.results.filt,
                    features = markers,
                    subset.nhoods = da.results$da.clusters %in% c(i),
                    assay = "logcounts",
                    scale_to_1 = TRUE,
                    cluster_features = TRUE,
                    show_rownames = FALSE,
                    group = c(i),
                    condition = condition,
                    alpha = fdr.cutoff
                  )
              
               PrintSaveAndSVG(p5, paste0("milo_DA_DE_heatmap_", i, "_", condition), plots.format)
              }
            }
          } else {
            print("no sig markers")
          }
        } else {
          print("skipped")
        }
      }
      
      # plots the Differential abundance -> uses UMAP reduction from seurat object
      # - nodes = neighborhoods
      #   - color = fold change
      #   - node layout -> based on index node in UMAP of single cells
      # - edges = number cells shared btw adjacent neighborhoods
      print("plotnhood")
      p3 <-
        plotNhoodGraphDA(sc.integrated.milo.traj,
                         da.results,
                         layout = "UMAP",
                         alpha = fdr.cutoff) +  # for my UMAP negative FC denoted CTRL and positive values denoted TREAT -> will have to look at your original UMAP for density of cells
      # also can look at the your design -> should be ordered correspondingly; https://github.com/MarioniLab/miloR/issues/81
        ggtitle(paste0("DA Analysis UMAP for ", condition))
      
      print("recalc da.results")
      PrintSaveAndSVG(p3, paste0("milo_DA_umap_", condition), plots.format)
      
      saveRDS(sc.integrated.milo.traj, paste0("sc_integrated_milo_traj_", condition, ".rds"))
    }
    print("setnull")
    se.integrated$da.clusters <- NULL
  }
