RunPsuper <- function(se.integrated, main_time, plots.format){
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
      
      beta_nzero 	= psuper.obj$beta_dt[ abs_beta > 0 ]
      n_nzero 	= nrow(beta_nzero)
      top_genes 	= as.character(beta_nzero[1:min(20, nrow(beta_nzero))]$symbol)
      
      if (nrow(beta_nzero) > 0) {
        p4 <- plot_identified_genes_over_psupertime(psuper.obj, palette = "Dark2") + ggtitle(paste0("Top 20 Genes over psupertime ", j))
  
        gene.comparison <- as.vector(levels(p4[["data"]][["symbol"]]))
        write.table(gene.comparison, paste0(data.dir, "psuper_top_20_genes_", j, ".txt"), quote = FALSE,row.names = F, sep = "\t", col.names = F)
      }
      
      ggSaveAndJPEG(p1,
                   paste0(tar.dir, "psuper_train_results_", j),
                   plots.format,
                   width = 8,
                   height = 6
                   )
      ggSaveAndJPEG(p2,
                   paste0(tar.dir, "psuper_density_pseudotime_", j),
                   plots.format,
                   width = 8,
                   height = 6
                   )
      ggSaveAndJPEG(p3,
                   paste0(tar.dir, "psuper_gene_coefficients_", j),
                   plots.format,
                   width = 8,
                   height = 6
                   )
      ggSaveAndJPEG(p4,
                   paste0(tar.dir, "psuper_top_20_genes_over_pseudotime_", j),
                   plots.format,
                   width = 8,
                   height = 6
                   )
    }
  } else {
    if(length(unique(se.integrated@meta.data[["group"]])) >= 2){
      density.summary <- list()
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
          
          ggSaveAndJPEG(p1,
                       paste0(tar.dir, "psuper_train_results_", j),
                       plots.format,
                       width = 8,
                       height = 6
          )
          ggSaveAndJPEG(p2,
                       paste0(tar.dir, "psuper_density_pseudotime_", j),
                       plots.format,
                       width = 8,
                       height = 6
          )
          ggSaveAndJPEG(p3,
                       paste0(tar.dir, "psuper_gene_coefficients_", j),
                       plots.format,
                       width = 8,
                       height = 6
          )
          ggSaveAndJPEG(p4,
                       paste0(tar.dir, "psuper_top_20_genes_over_pseudotime_", j),
                       plots.format,
                       width = 8,
                       height = 6
          )
          
          psuper.df <- psuper.obj[["proj_dt"]]
          psuper.df$group <- i
          psuper.df$cluster <- j
          density.summary <- append(density.summary, list(psuper.df))
          
          for (k in unique(se.integrated@meta.data[["group"]])){
            if (i != k){
              se.integrated.filt <- subset(se.integrated, group == i & psupertime.labels == j)
              y <- se.integrated.filt@meta.data[["time"]]
              levels(y) <- mixedsort(unique(y))
              fac_y <- factor(y, levels = mixedsort(unique(y)))
              psuper.obj2  = psupertime(se.integrated.filt@assays$RNA$counts, fac_y, sel_genes='list', gene_list = as.vector(gene.comparison), penalization='best') # heurstics
              
              p3 <- plot_identified_gene_coefficients(psuper.obj2) + ggtitle(paste0(i, ": ", "Top Gene Coefficients ", j, " in ", k))
              p4 <- plot_identified_genes_over_psupertime(psuper.obj2, palette = "Dark2") + ggtitle(paste0(i, ": ", "Top 20 Genes over psupertime ", j, " in ", k))
              
              ggSaveAndJPEG(p3,
                           paste0(tar.dir, "psuper_gene_coefficients_", i, "_", j, "_", k, "_genes"),
                           plots.format,
                           width = 8,
                           height = 6
                           )
              ggSaveAndJPEG(p4,
                           paste0(tar.dir, "psuper_top_20_genes_over_pseudotime_", i, "_", j, "_", k, "_genes"),
                           plots.format,
                           width = 8,
                           height = 6
                           )
            }
          }
        }
      }
      density.summary.bind <- bind_rows(density.summary)
      for (i in unique(density.summary.bind$cluster)) {
        density.summary.bind.filt <- subset(density.summary.bind, cluster == i)
        stat.test <- density.summary.bind.filt %>%
          group_by(label_input) %>%
          t_test(psuper ~ group) %>%
          adjust_pvalue(method = "BH") %>%
          add_significance()
        
        # ggplot(density.summary.bind.filt, aes(x=label_input, y=psuper , fill=group)) +
        #   geom_boxplot() + stat_pvalue_manual(stat.test, y.position = 1, label = "p.adj.signif", remove.bracket = T)
        
        for (j in 1:length(rownames(stat.test))){
          stat.test[j, 3:4] <- stat.test[j, 1]
        }
        
        p <- ggboxplot(density.summary.bind.filt, x = "label_input", y = "psuper",
                       color = "group") +
          ggtitle(paste0(i, ": Psupertime Distribution")) + 
          theme(legend.position = "right",
                axis.title.x=element_blank(),
                axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          ylab("psupertime")
        
        p <- p + stat_pvalue_manual(stat.test, y.position = max(density.summary.bind.filt$psuper) + 0.1, label = "p.adj.signif", remove.bracket = T)
        ggSaveAndJPEG(p,
                     paste0("psupertime_plots/", "psuper_boxplot_compare_dist_", i),
                     plots.format,
                     width = 8,
                     height = 6
                     )
        
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
        
        ggSaveAndJPEG(p1,
                     paste0(tar.dir, "psuper_train_results_", j),
                     plots.format,
                     width = 8,
                     height = 6
        )
        ggSaveAndJPEG(p2,
                     paste0(tar.dir, "psuper_density_pseudotime_", j),
                     plots.format,
                     width = 8,
                     height = 6
        )
        ggSaveAndJPEG(p3,
                     paste0(tar.dir, "psuper_gene_coefficients_", j),
                     plots.format,
                     width = 8,
                     height = 6
        )
        ggSaveAndJPEG(p4,
                     paste0(tar.dir, "psuper_top_20_genes_over_pseudotime_", j),
                     plots.format,
                     width = 8,
                     height = 6
        )
      }
    }
  }
}