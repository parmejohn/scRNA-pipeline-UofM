set.seed(333)

#setwd("/research/2024_scrnaseq_pipeline/organoid_work/cellchat_testing_mult_cond")

CellChatAnalysis <- function(se.integrated, species){
  if(any(colnames(se.integrated[[]]) == "celltype")){
    Idents(se.integrated) <- se.integrated$celltype
  } else if (any(colnames(se.integrated[[]]) == "wsnn_res.1")) {
    se.integrated$cluster_wsnn<- paste0("cluster_", se.integrated$wsnn_res.1)
    Idents(se.integrated) <- se.integrated$cluster_wsnn
  } else {
    se.integrated$cluster_seurat <- paste0("cluster_", se.integrated$seurat_clusters)
    Idents(se.integrated) <- se.integrated$cluster_seurat
  }
  
  
  ##### Setting up database #####
  CellChatDB <- NA
  if (species == "Homo sapiens"){
    CellChatDB <- CellChatDB.human
  } else if (species == "Mus musculus"){
    CellChatDB <- CellChatDB.mouse
  }
  CellChatDB.use <- CellChatDB # using all signals = metabolic and synaptic as well
  
  ##### Creating and processing cellchat data (the receptor-ligand information) #####
  group.pairs <- NA
  list.comparisons <- list("group")
  if (length(se.integrated@misc) != 0 ){
    for (j in se.integrated@misc$co.conditions){
      group.co <- paste0("group_", j)
      se.integrated <- AddMetaData(se.integrated, paste(se.integrated$group, se.integrated@meta.data[[j]], sep = "_"), group.co)
      list.comparisons <- append(list.comparisons, group.co)
    }
  }
  
  dir.create("cellchat_data")
  dir.create("cellchat_plots")
  
  for (k in 1:length(list.comparisons)){
    #Idents(se.integrated) <- k
    grouping <- list.comparisons[[k]]
    if (any(duplicated(as.data.frame(se.integrated@meta.data[[grouping]])))){ # check if you have at least more than 1 sample for the comparison
      
      # create pairs for each possible combination
      group.pairs <- as.data.frame(combn(unique(se.integrated@meta.data[[grouping]]), 2))
      group.pairs <- sapply(group.pairs, function(x) as.character(x), simplify = FALSE)
      #print(group.pairs)
      #if (sub("^[^_]*_", "", k) == "time"){
      if (k >= 2){
        group.pairs.filt <- list()
        for (z in group.pairs){
          time.match <- sub("^[^_]*_", "", z)
          if (time.match[1] == time.match[2]) {
            cond.mismatch <- sub("_.*", "", z)
            if (cond.mismatch[1] != cond.mismatch[2]){
              group.pairs.filt <- append(group.pairs.filt, list(z))
            }
          }
        }
        group.pairs <- group.pairs.filt
      }
    }
    for (z in group.pairs){
      data.dir <- paste0("cellchat_data/", z[1], "_vs_", z[2], "/")
      dir.create(data.dir)
      
      plots.dir <- paste0("cellchat_plots/", z[1], "_vs_", z[2], "/")
      dir.create(plots.dir)
      
      object.list <- list()
      for (i in z){
        # set up dataframe
        print(i)
        expr <- FetchData(se.integrated, vars = grouping)
        se.integrated.group <- se.integrated[, which(expr == i)]
        
        #se.integrated.group <- subset(se.integrated, group == i)
        
        data.input <- se.integrated.group[["RNA"]]$data
        labels <- Idents(se.integrated.group)
        samples <- se.integrated.group$sample
        
        meta <- data.frame(labels = labels, samples = samples, row.names = names(labels)) # create a dataframe of the cell labels
        
        cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
        cellchat@DB <- CellChatDB.use
        
        cellchat <- subsetData(cellchat)
        cellchat <- identifyOverExpressedGenes(cellchat)
        cellchat <- identifyOverExpressedInteractions(cellchat)
        
        cellchat <- computeCommunProb(cellchat, type = "triMean")
        cellchat <- filterCommunication(cellchat, min.cells = 10)
        cellchat <- computeCommunProbPathway(cellchat)
        cellchat <- aggregateNet(cellchat)
        cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
        
        object.list <- c(object.list, cellchat)
      }
      
      names(object.list) <-  unique(z)
      saveRDS(object.list, paste0(data.dir, "cellchat_object_list.rds"))
      
      cellchat.merged <- mergeCellChat(object.list, add.names = names(object.list))
      saveRDS(cellchat.merged, paste0(data.dir, "cellchat_merged.rds"))
      #cellchat.merged <- readRDS("~/cellchat_merged.rds")
      
      ## compares the number of interactions and the strength from given conditions
      gg1 <- compareInteractions(cellchat.merged, show.legend = F, group = c(1,2))
      gg2 <- compareInteractions(cellchat.merged, show.legend = F, group = c(1,2), measure = "weight")
      
      p1 <- gg1 + gg2
      ggsave(paste0(plots.dir, "cellchat_interaction_summary_bar.pdf"), p1, width = 8, height = 6)
      
      # red = increased signaling in the second dataset compared to the first
      pdf(paste0(plots.dir, "cellchat_differential_interaction_circle.pdf"), width = 10, height = 6)
      par(mfrow = c(1,2), xpd=TRUE)
      tmp1 <- netVisual_diffInteraction(cellchat.merged, weight.scale = T, margin = 0.3)
      tmp2 <- netVisual_diffInteraction(cellchat.merged, weight.scale = T, measure = "weight", margin = 0.3)
      graphics.off()
      
      # top colored bar plot -> sum of each column for incoming signaling
      # right colored bar plot = sum of each row -> outgoing signal
      gg1 <- netVisual_heatmap(cellchat.merged)
      #> Do heatmap based on a merged object
      gg2 <- netVisual_heatmap(cellchat.merged, measure = "weight")
      #> Do heatmap based on a merged object
      p2 <- draw(gg1 + gg2, ht_gap = unit(0.5, "cm"))
      
      # pdf(paste0(plots.dir, "cellchat_differential_interaction_heatmap.pdf"), width = 8, height = 6)
      # p2
      # graphics.off()
      PrintSave(p2, paste0(plots.dir, "cellchat_differential_interaction_heatmap.pdf"), w = 8, h = 6)
      
      # Compare num of interactions or interaction strength across cell populations
      weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
      
      pdf(paste0(plots.dir, "cellchat_num_interactions_circle.pdf"), width = 10, height = 6)
      par(mfrow = c(1,2), xpd=TRUE)
      for (i in 1:length(object.list)) {
        netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], 
                         edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
      }
      graphics.off()
      
      
      # find cell populations with sig chgs in sending/receiving signals
      num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
      weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
      gg <- list()
      for (i in 1:length(object.list)) {
        gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
      }
      #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
      #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
      p3 <- patchwork::wrap_plots(plots = gg)
      ggsave(paste0(plots.dir, "cellchat_population_send_receive.pdf"), p3, width = 8, height = 6)
      
      ##### Find altered signalling with distinct interaction strength (also called information flow) #####
      # information flow = sum of communication prob among all pairs of cell groups
      # sig pathways -> based on diff in overall info flow
      # paired wilcox test is used to see if there is a sig diff of signalling info
      gg1 <- rankNet(cellchat.merged, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
      gg2 <- rankNet(cellchat.merged, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
      
      p4 <- gg1 + gg2
      ggsave(paste0(plots.dir, "cellchat_information_flow_compare.pdf"), p4, width = 8, height = 6)
      
      ##### get the top signalling pathways and plot out cirle plots #####
      top.paths <- subset(gg1[["data"]], pvalues <= 0.05)
      top.paths <- unique(top.paths$name)
      
      cellchat.merged@meta$datasets = factor(cellchat.merged@meta$datasets, levels = c(z[1], z[2])) # set factor level
      
      dir.create(paste0(plots.dir, "signaling_pathways"))
      for (i in top.paths){
        pathways.show <- c(i)
        dir.create(paste0(plots.dir, "signaling_pathways/", i))
        tar.dir <- paste0(plots.dir, "signaling_pathways/", i)
        if (i %in% object.list[[z[1]]]@netP[["pathways"]] & i %in% object.list[[z[2]]]@netP[["pathways"]]) {
          weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
          # par(mfrow = c(1,2), xpd=TRUE)
          # for (j in 1:length(object.list)) {
          #   netVisual_aggregate(object.list[[j]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[j]))
          # }
          pdf(file = paste0(tar.dir, "/cellchat_", i, "_signal_path.pdf"), width = 10, height = 10)
          par(mfrow = c(1,2), xpd=TRUE)
          for (j in 1:length(object.list)) {
            netVisual_aggregate(object.list[[j]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[j]))
          }
          graphics.off()
        }
        p5 <- plotGeneExpression(cellchat.merged, signaling = i, split.by = "datasets", colors.ggplot = T, type = "violin") +
          plot_annotation(paste0(i, " pathway gene expression"),theme=theme(plot.title=element_text(hjust=0.5)))
        ggsave(paste0(tar.dir, "/cellchat_", i, "_expression.pdf"), p5, width = 8, height = 8)
      }
      
      ##### Outgoing/incoming signalling patterns for each cell pop #####
      i = 1
      # combining all the identified signaling pathways from different datasets 
      pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
      ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
      ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
      p6 <- draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
      PrintSave(p6, "cellchat_compare_outgoing_signal_heatmap.pdf", plots.dir, w = 8, h = 8)
      
      ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
      ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
      p7 <- draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
      PrintSave(p7, "cellchat_compare_incoming_signal_heatmap.pdf", plots.dir, w = 8, h = 8)
      
      ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
      ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
      p8 <- draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
      PrintSave(p8, "cellchat_compare_all_signal_heatmap.pdf", plots.dir, w = 8, h = 8)
      
      ##### identifying LR pairs #####
      # this requires specific cell types -> get the number of incoming/outgoing interaction strength and get the top 3 cell types from there -> too many for 1 graph still
      # use a for loop and print all results
      # make a plot of each celltype
      # pval = significance levels of communication, not btw the 2 conditions
      # communication prob = interaction strength increasing/decreasing
      ## Graph is filtered by the missing prescence of the communication prob in the other condition, or if there is an abs(log2FC) >= 0.6 -> ~50% increase in prob
      ## Sometimes the dataset cannot be filtered (NA source targets after filtering), so print the unfiltered option in these situations
      dir.create(paste0(plots.dir, "commun_prob"))
      for (i in 1:length(levels(se.integrated.group$ident))){
        
        se.integrated.check <- subset(se.integrated.group, ident == levels(se.integrated.group$ident)[i])
        min.cells <- ncol(se.integrated.check@assays[["RNA"]])
        
        if (min.cells > 10){
          gg1 <- netVisual_bubble(cellchat.merged, sources.use = levels(se.integrated.group$ident)[i], 
                                  comparison = c(1, 2), max.dataset = 2, title.name = "", angle.x = 45, remove.isolate = T, return.data = T)
          
          gg1 <- filterLRPairsBubble(cellchat.merged, gg1, i, max.dataset = 2, title.name = paste0("Increased signaling in ", z[2]), ident.1 = z[1], ident.2 = z[2])
          
          gg2 <- netVisual_bubble(cellchat.merged, sources.use = levels(se.integrated.group$ident)[i],
                                  comparison = c(1, 2), max.dataset = 1, title.name = "", angle.x = 45, remove.isolate = T, return.data = T)
          
          gg2 <- filterLRPairsBubble(cellchat.merged, gg2, i, max.dataset = 1, title.name = paste0("Decreased signaling in ", z[2]), ident.1 = z[1], ident.2 = z[2])
          
          p9 <- gg1 + gg2
          ggsave(paste0(plots.dir, "commun_prob/cellchat_", levels(se.integrated.group$ident)[i], "_expression.pdf"), p9, width = 12, height = 8)
        }
      }
    }
  }


}

filterLRPairsBubble <- function(cellchat.merged, gg1, i, max.dataset, title.name, ident.1, ident.2){
  #gg1 <- netVisual_bubble(cellchat.merged, sources.use = i,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in WT", angle.x = 45, remove.isolate = T, return.data = T)
  
  ## row numbers that are uniquely WT or MUT
  unique_row_numbers <- gg1[["gg.obj"]][["data"]] %>%
    mutate(row_number = row_number()) %>%  # Add a column with row numbers
    group_by(interaction_name, group.names) %>%  # Group by columns that identify a pair
    filter(n_distinct(dataset) == 1) %>%  # Keep groups that have only one unique value in V16
    ungroup() %>%
    pull(row_number)  # Extract the row numbers
  
  ## find those with a log2FC >= 1 in probability if the MUT and WT pair exists
  log2fc_results_rows <- gg1[["gg.obj"]][["data"]] %>%
    mutate(row_number = row_number()) %>%  # Add a column with row numbers
    group_by(interaction_name, group.names) %>%
    filter(n() == 2 & all(c(ident.2, ident.1) %in% dataset)) %>%  # Ensure both WT and MUT are present
    summarize(
      log2FC = log2(prob[dataset == ident.1] / prob[dataset == ident.2]),
      ident_1_row = row_number[dataset == ident.1],
      ident_2_row = row_number[dataset == ident.2]
    ) %>%
    ungroup() %>%
    filter(abs(log2FC) >= 0.58) %>%  # Filter rows where log2FC is greater than 0.58
    select(ident_1_row, ident_2_row) %>%  # Select the row numbers
    unlist()
  
  to_filter <- c(unique_row_numbers, log2fc_results_rows)
  
  tmp <-  gg1[["gg.obj"]][["data"]]
  gg1[["gg.obj"]][["data"]] <- gg1[["gg.obj"]][["data"]][to_filter,]
  
  if (length(unique(gg1[["gg.obj"]][["data"]]$interaction_name)) != 0){
    interaction_name.filtered <- as.data.frame(unique(gg1[["gg.obj"]][["data"]]$interaction_name))
    colnames(interaction_name.filtered)[1] <- "interaction_name"
    if (length(which(is.na(as.character(gg1[["gg.obj"]][["data"]]$source.target)))) == 0 & length(gg1[["gg.obj"]][["data"]]$source.target) >= 4){
      p1 <- netVisual_bubble(cellchat.merged, sources.use = i, pairLR.use = interaction_name.filtered, comparison = c(1, 2), 
                             max.dataset = max.dataset, title.name = title.name, angle.x = 45, remove.isolate = T)
      return(p1)
    } else {
      #for (j in which(is.na(as.character(gg1[["gg.obj"]][["data"]]$source.target)))){
      #levels(gg1[["gg.obj"]][["data"]]$source.target) <- factor(c(levels(gg1[["gg.obj"]][["data"]]$source.target), paste0(gg1[["gg.obj"]][["data"]][i,14], " ", "(", gg1[["gg.obj"]][["data"]][j,15], ")")))
      #gg1[["gg.obj"]][["data"]][j, 12] <- paste0(gg1[["gg.obj"]][["data"]][j,14], " ", "(", gg1[["gg.obj"]][["data"]][j,15], ")")
      p1 <- netVisual_bubble(cellchat.merged, sources.use = i, comparison = c(1, 2), max.dataset = max.dataset, title.name = paste0(title.name, ": Unfiltered"), angle.x = 45, remove.isolate = T)
      return(p1)
    }
  } 
}
