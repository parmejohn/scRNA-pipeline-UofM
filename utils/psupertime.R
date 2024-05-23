RunPsuper <- function(se.integrated, main_time){
  dir.create("psupertime_plots")
  dir.create("psupertime_data")
  se.integrated$psupertime.labels <- Idents(se.integrated)
  
  ## if time == main condition then take the grouping variable instead
  if (main_time == "yes"){
    for (j in levels(se.integrated$psupertime.labels)){
      dir.create(paste0("psupertime_plots/", j, "/"))
      tar.dir <- paste0("psupertime_plots/", j, "/")
      dir.create(paste0("psupertime_data/", j, "/"))
      data.dir <- paste0("psupertime_data/", j, "/")
      
      se.integrated.filt <- subset(se.integrated, psupertime.labels == j)
      y <- se.integrated.filt@meta.data[["group"]]
      levels(y) <- mixedsort(unique(y))
      fac_y <- factor(y, levels = mixedsort(unique(y)))
      psuper.obj  = psupertime(se.integrated.filt@assays$RNA$counts, fac_y, sel_genes='hvg', penalization='best') # heurstics

      p1 <- plot_train_results(psuper.obj) + ggtitle(paste0("Train Results ", j))
      p2 <- plot_labels_over_psupertime(psuper.obj) + ggtitle(paste0("Labels over Time ", j))
      p3 <- plot_identified_gene_coefficients(psuper.obj) + ggtitle(paste0("Top Gene Coefficients ", j))
      p4 <- plot_identified_genes_over_psupertime(psuper.obj, palette = "Dark2") + ggtitle(paste0("Top 20 Genes over psupertime ", j))

      gene.comparison <- as.vector(levels(p4[["data"]][["symbol"]]))
      write.table(gene.comparison, paste0(data.dir, "psuper_top_20_genes_", j, ".txt"), quote = FALSE,row.names = F, sep = "\t", col.names = F)
      
      ggsave(paste0(tar.dir, "psuper_train_results_", j, ".pdf"), plot = p1, width = 8, height = 6)
      ggsave(paste0(tar.dir, "psuper_density_pseudotime_", j, ".pdf"), plot = p2, width = 8, height = 6)
      ggsave(paste0(tar.dir, "psuper_gene_coefficients_", j, ".pdf"), plot = p3, width = 8, height = 6)
      ggsave(paste0(tar.dir, "psuper_top_20_genes_over_pseudotime_", j, ".pdf"), plot = p4, width = 8, height = 6)
    }
  } else {
    if(length(unique(se.integrated@meta.data[["group"]])) >= 2){
      for (i in unique(se.integrated@meta.data[["group"]])){
        dir.create(paste0("psupertime_plots/", i))
        dir.create(paste0("psupertime_data/", i))
        #gene.comparison <- NA # gene list will change for each cell cluster though
        for (j in levels(se.integrated$psupertime.labels)){
          dir.create(paste0("psupertime_plots/", i, "/", j, "/"))
          tar.dir <- paste0("psupertime_plots/", i, "/", j, "/")
          dir.create(paste0("psupertime_data/",  i, "/", j, "/"))
          data.dir <- paste0("psupertime_data/",  i, "/", j, "/")
          
          se.integrated.filt <- subset(se.integrated, group == i & psupertime.labels == j)
          y <- se.integrated.filt@meta.data[["time"]]
          levels(y) <- mixedsort(unique(y))
          fac_y <- factor(y, levels = mixedsort(unique(y)))
          psuper.obj  = psupertime(se.integrated.filt@assays$RNA$counts, fac_y, sel_genes='hvg', penalization='best') # heurstics

          p1 <- plot_train_results(psuper.obj) + ggtitle(paste0(i, ": ", "Train Results ", j))
          p2 <- plot_labels_over_psupertime(psuper.obj) + ggtitle(paste0(i, ": ", "Labels over Time ", j))
          p3 <- plot_identified_gene_coefficients(psuper.obj) + ggtitle(paste0(i, ": ", "Top Gene Coefficients ", j))
          p4 <- plot_identified_genes_over_psupertime(psuper.obj, palette = "Dark2") + ggtitle(paste0(i, ": ", "Top 20 Genes over psupertime ", j))
          
          gene.comparison <- as.vector(levels(p4[["data"]][["symbol"]]))
          write.table(gene.comparison, paste0(data.dir, "psuper_top_20_genes_", i, "_", j, ".txt"), quote = FALSE,row.names = F, sep = "\t", col.names = F)
          
          ggsave(paste0(tar.dir, "psuper_train_results_", i, "_", j, ".pdf"), plot = p1, width = 8, height = 6)
          ggsave(paste0(tar.dir, "psuper_density_pseudotime_", i, "_", j, ".pdf"), plot = p2, width = 8, height = 6)
          ggsave(paste0(tar.dir, "psuper_gene_coefficients_", i, "_", j, ".pdf"), plot = p3, width = 8, height = 6)
          ggsave(paste0(tar.dir, "psuper_top_20_genes_over_pseudotime_", i, "_", j, ".pdf"), plot = p4, width = 8, height = 6)
          
          for (k in unique(se.integrated@meta.data[["group"]])){
            if (i != k){
              se.integrated.filt <- subset(se.integrated, group == i & psupertime.labels == j)
              y <- se.integrated.filt@meta.data[["time"]]
              levels(y) <- mixedsort(unique(y))
              fac_y <- factor(y, levels = mixedsort(unique(y)))
              psuper.obj  = psupertime(se.integrated.filt@assays$RNA$counts, fac_y, sel_genes='list', gene_list = as.vector(gene.comparison), penalization='best') # heurstics
              
              p3 <- plot_identified_gene_coefficients(psuper.obj) + ggtitle(paste0(i, ": ", "Top Gene Coefficients ", j, " in ", k))
              p4 <- plot_identified_genes_over_psupertime(psuper.obj, palette = "Dark2") + ggtitle(paste0(i, ": ", "Top 20 Genes over psupertime ", j, " in ", k))
              
              # ggsave(paste0(tar.dir, "psuper_train_results_", i, "_", j, "_", k, "_genes", ".pdf"), plot = p1, width = 8, height = 6)
              # ggsave(paste0(tar.dir, "psuper_density_pseudotime_", i, "_", j, "_", k, "_genes", ".pdf"), plot = p2, width = 8, height = 6)
              ggsave(paste0(tar.dir, "psuper_gene_coefficients_", i, "_", j, "_", k, "_genes", ".pdf"), plot = p3, width = 8, height = 6)
              ggsave(paste0(tar.dir, "psuper_top_20_genes_over_pseudotime_", i, "_", j, "_", k, "_genes", ".pdf"), plot = p4, width = 8, height = 6)    
            }
          }
        }
      }
    } else {
      for (j in levels(se.integrated$psupertime.labels)){
        dir.create(paste0("psupertime_plots/", j, "/"))
        tar.dir <- paste0("psupertime_plots/", j, "/")
        dir.create(paste0("psupertime_data/", j, "/"))
        data.dir <- paste0("psupertime_data/", j, "/")
        
        se.integrated.filt <- subset(se.integrated, psupertime.labels == j)
        y <- se.integrated.filt@meta.data[["time"]]
        levels(y) <- mixedsort(unique(y))
        fac_y <- factor(y, levels = mixedsort(unique(y)))
        psuper.obj  = psupertime(se.integrated.filt@assays$RNA$counts, fac_y, sel_genes='hvg', penalization='best') # heurstics
        
        p1 <- plot_train_results(psuper.obj) + ggtitle(paste0("Train Results ", j))
        p2 <- plot_labels_over_psupertime(psuper.obj) + ggtitle(paste0("Labels over Time ", j))
        p3 <- plot_identified_gene_coefficients(psuper.obj) + ggtitle(paste0("Top Gene Coefficients ", j))
        p4 <- plot_identified_genes_over_psupertime(psuper.obj, palette = "Dark2") + ggtitle(paste0("Top 20 Genes over psupertime ", j))
        
        gene.comparison <- as.vector(levels(p4[["data"]][["symbol"]]))
        write.table(gene.comparison, paste0(data.dir, "psuper_top_20_genes_", j, ".txt"), quote = FALSE,row.names = F, sep = "\t", col.names = F)
        
        ggsave(paste0(tar.dir, "psuper_train_results_", j, ".pdf"), plot = p1, width = 8, height = 6)
        ggsave(paste0(tar.dir, "psuper_desnsity_pseudotime_", j, ".pdf"), plot = p2, width = 8, height = 6)
        ggsave(paste0(tar.dir, "psuper_gene_coefficients_", j, ".pdf"), plot = p3, width = 8, height = 6)
        ggsave(paste0(tar.dir, "psuper_top_20_genes_over_pseudotime_", j, ".pdf"), plot = p4, width = 8, height = 6)
      }
    }
  }
  
}