set.seed(333)

jpeg <- function(filename, width, height, dpi = 300) {
	  # Convert inches to pixels
	  pixel_width <- width * dpi
  pixel_height <- height * dpi
    
    # Call the original jpeg function
    grDevices::jpeg(filename, width = pixel_width, height = pixel_height, res = dpi)
}

PrintSave <- function(plot, title, plots.format = "pdf", width=8, height=6){
  title.with.ext <- paste0(title, ".", plots.format)
  plot_function <- get(plots.format)
  
  plot_function(title.with.ext, width = width, height = width)
  print(plot)
  graphics.off()
}

PrintSaveAndJPEG <- function(plot, title, plots.format, width=8, height=6){
  PrintSave(plot, title, plots.format, width = width, height = height)
  PrintSave(plot, title, "jpeg", width = width, height = height)
}

ggSaveAndJPEG <- function(plot, title, plots.format, width=8, height=8){
  ggsave(paste0(title, ".", plots.format), plot, width = width, height = height)
  ggsave(paste0(title, ".jpeg"), plot, width = width, height = height)
  
}

ListAllPossibleComparisons <- function(se.integrated, seurat.subset) {
  list.comparisons <- list("group")
  if (length(se.integrated@misc) != 0){
    for (j in se.integrated@misc$co.conditions){
      group.co <- paste0("group_", j)
      seurat.subset <- AddMetaData(seurat.subset, 
                                             paste(seurat.subset$group, 
                                                   seurat.subset@meta.data[[j]], sep = "_"), 
                                             group.co
                                             )
      list.comparisons <- append(list.comparisons, group.co)
    }
  }
  return(list(subset = seurat.subset, comparisons = list.comparisons))
}

MatchCovariantGroupPairs <- function(seurat.subset, grouping, not.main.group) {
  # create pairs for each possible combination
  group.pairs <- as.data.frame(combn(unique(seurat.subset@meta.data[[grouping]]), 2))
  group.pairs <- sapply(group.pairs, function(x) as.character(x), simplify = FALSE)
  
  if (not.main.group >= 2){
    group.pairs.filt <- list()
    for (z in group.pairs){
      time.match <- sub("^[^_]*_", "", z) # not just time, this is for any co-variant
      if (time.match[1] == time.match[2]) {
        cond.mismatch <- sub("_.*", "", z)
        if (cond.mismatch[1] != cond.mismatch[2]){
          group.pairs.filt <- append(group.pairs.filt, list(z))
        }
      }
    }
    group.pairs <- group.pairs.filt
    print(group.pairs)
    
  }
  return(group.pairs)
}

Correcth5SeuratFile <- function(file){
  f <- H5File$new(file, "r+")
  groups <- f$ls(recursive = TRUE)
  
  for (name in groups$name[grepl("categories", groups$name)]) {
    names <- strsplit(name, "/")[[1]]
    names <- c(names[1:length(names) - 1], "levels")
    new_name <- paste(names, collapse = "/")
    f[[new_name]] <- f[[name]]
  }
  
  for (name in groups$name[grepl("codes", groups$name)]) {
    names <- strsplit(name, "/")[[1]]
    names <- c(names[1:length(names) - 1], "values")
    new_name <- paste(names, collapse = "/")
    f[[new_name]] <- f[[name]]
    grp <- f[[new_name]]
    grp$write(args = list(1:grp$dims), value = grp$read() + 1)
  }
  
  f$close_all()
}
