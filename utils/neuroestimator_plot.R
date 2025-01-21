NeuroestimatorPlot <- function(se.integrated, neuroestimator.results, plots.format){
	
	se.integrated <- readRDS(se.integrated)
  neuroestimator.results <- read.table(neuroestimator.results)
  metadata <- se.integrated@meta.data

  if ("celltype" %in% names(metadata)){
    anno <- "celltype"
  } else {
    anno <- "seurat_clusters"
  }
  
  # Subset metadata to include only the relevant column(s) and matching rows
  if (length(se.integrated@misc[["co.conditions"]]) > 0) {
    metadata_subset <- metadata[, c("group", anno, se.integrated@misc[["co.conditions"]]), drop = FALSE]  # Replace with your column names
  } else {
    metadata_subset <- metadata[, c("group"), anno, drop = FALSE]  # Replace with your column names
  }
  metadata_subset <- metadata_subset[rownames(neuroestimator.results), , drop = FALSE]  # Match barcodes
  
  # Combine predicted_activity with metadata
  neuroestimator.results.meta <- cbind(neuroestimator.results, metadata_subset)
  
  if (length(se.integrated@misc[["co.conditions"]]) > 0) {
    for (i in se.integrated@misc[["co.conditions"]]){
      co.variant <- select(neuroestimator.results.meta, 'predicted_activity', 'group', anno, i)
      extra_column <- i
      stat.test <- KSTestOverGroups(anno, extra_column, co.variant)
      PlotNeuroestimatorResults(stat.test, co.variant, extra_column, anno, plots.format)
      graphics.off()
    }
  } else {
    co.variant <- neuroestimator.results.meta
    extra_column <- ""
    stat.test <- KSTestOverGroups(anno, extra_column, co.variant)
    PlotNeuroestimatorResults(stat.test, co.variant, extra_column, anno, plots.format)
    graphics.off()
  }
}

KSTestOverGroups <- function(anno, extra_column, co.variant) {
  group_vars <- c(anno)
  if (extra_column %in% colnames(co.variant)) {
    group_vars <- c(group_vars, extra_column)
  }
  
  stat.test <- co.variant %>%
    group_by(across(all_of(group_vars))) %>%  # Dynamically group by both columns if present
    summarise(
      ks_result = list(ks.test(predicted_activity[group == unique(group)[1]], 
                               predicted_activity[group == unique(group)[2]])),
      .groups = "drop"
    ) %>%
    mutate(
      p_value = sapply(ks_result, function(x) x$p.value),
      statistic = sapply(ks_result, function(x) x$statistic)
    ) %>%
    select(-ks_result) %>%
    mutate(
      p_value_adjusted = p.adjust(p_value, method = "BH"),
      significance = case_when(
        p_value_adjusted < 0.001 ~ "***",
        p_value_adjusted < 0.01 ~ "**",
        p_value_adjusted < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
}

PlotNeuroestimatorResults <- function(stat.test, co.variant, extra_column, anno, plots.format) {
  ident.1 <- unique(co.variant$group)[1]
  ident.2 <- unique(co.variant$group)[2]
  
  if (extra_column != ""){
    for (i in unique(stat.test[,2])[[1]]) {
      p <- co.variant %>% 
        filter(!!sym(extra_column) == i) %>%
        ggboxplot(x = eval(anno), 
                  y = "predicted_activity",
                  color = "group") +
        ylim(0,1) +
        ggtitle(paste0("Predicted neuronal activity for ", ident.1, " vs ", ident.2, ": ", i)) + 
        theme(legend.position = "right",
              axis.title.x=element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        ylab("NEUROeSTIMator Score") 
      
      stat.test.filtered <- filter(stat.test %>%
                                     filter(!!sym(extra_column) == i))
      
      stat.test.filtered$group1 <- NA
      stat.test.filtered$group2 <- NA
      
      for (j in 1:length(rownames(stat.test.filtered))){
        stat.test.filtered[j, 7] <- stat.test.filtered[j, 1]
        stat.test.filtered[j, 8] <- stat.test.filtered[j, 1]
      } 
      
      p <- p + stat_pvalue_manual(stat.test.filtered, y.position = 1, label = "significance", remove.bracket = T)
      
      ggSaveAndJPEG(plot = p, paste0("neuroestimator_boxplot_", ident.1, "_vs_", ident.2, "_", i), plots.format, width = 8, height = 8)
    }
  } else {
    p <- co.variant %>% 
      ggboxplot(x = eval(anno), 
                y = "predicted_activity",
                color = "group") +
      ylim(0,1) +
      ggtitle(paste0("Predicted neuronal activity for ", ident.1, " vs ", ident.2)) + 
      theme(legend.position = "right",
            axis.title.x=element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ylab("NEUROeSTIMator Score") 
    
    stat.test$group1 <- NA
    stat.test$group2 <- NA
    
    for (j in 1:length(rownames(stat.test))){
      stat.test[j, 6] <- stat.test[j, 1]
      stat.test[j, 7] <- stat.test[j, 1]
    } 
    
    p <- p + stat_pvalue_manual(stat.test, y.position = 1, label = "significance", remove.bracket = T)
    
    ggSaveAndJPEG(plot = p, paste0("neuroestimator_boxplot_", ident.1, "_vs_", ident.2), plots.format,width = 8, height = 8)
  }
}
