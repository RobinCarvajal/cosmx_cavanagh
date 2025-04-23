# load brain/celltypes_bulk objects
celltypes_bulk <- readRDS(file.path(data_dir,"broad_single_cell",
                                  "broad.labels_pseudo_bulk_condition_astro+immune.RDS"))

# creating a list of em tables for all celltypes -------------------------------
em_list <- list()

for (celltype in names(celltypes_bulk)){
  
  cluster_df <- as.data.frame(celltypes_bulk[[celltype]])
  
  cluster_df$ID <- rownames(cluster_df)
  
  cluster_df <- cluster_df[, c("ID", setdiff(names(cluster_df), "ID"))]
  
  # add to the em_list
  em_list[[celltype]] <- cluster_df
  
}



# creating a ss list for all cell types -------------------------------------------------
ss_list <- list()

for (celltype in names(celltypes_bulk)){
  
  # get the celltype
  counts <- celltypes_bulk[[celltype]]
  
  # generate sample level metadata
  ss <- tibble(sample = colnames(counts))
  
  # position
  segments <- c(2)
  replacement <- paste0("\\", segments[1])
  ss$position <- gsub('^([^_]+)_([^_]+)+_([^_]+)_([^_]+)$',
                       replacement, 
                       colnames(counts))
  # condition
  segments <- c(4)
  replacement <- paste0("\\", segments[1])
  ss$condition <- gsub('^([^_]+)_([^_]+)+_([^_]+)_([^_]+)$',
                          replacement, 
                          colnames(counts))
  
  
  
  # print to check
  print(ss)
  
  ss_list[[celltype]] <- ss
  
}

# ALL comparisons  -------------------------------------------------------

# NOTE: I do all comparisons in a loop to avoid repeating code

# in a loop
for (celltype in names(celltypes_bulk)) {
  
  # current celltype
  print(celltype)
  
  # Wrap in tryCatch to handle errors
  tryCatch({
    # base path
    base_path <- file.path(data_dir,"broad_single_cell","differential_expression", "position+condition")
    
    # out path 
    celltype_path <- file.path(base_path, celltype)
    # make dir
    if (!dir.exists(celltype_path)) {dir.create(celltype_path, recursive = TRUE)}
    
    # em 
    em <- em_list[[celltype]]
    # save the em 
    write.table(em, 
                file = file.path(celltype_path, paste0("em_", celltype, ".tsv")), 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
    # ss table
    ss <- ss_list[[celltype]]
    # save the ss
    write.table(ss,
                file = file.path(celltype_path, paste0("ss_", celltype, ".tsv")),
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
    
    # Control vs Aldara
    # GetDE_fn(em_df = em,
    #          var = celltype,
    #          suffix = "C_vs_A", 
    #          sample_sheet = ss,
    #          reference = "C", 
    #          contrast = "A", 
    #          p_threshold = 0.5, 
    #          fold_threshold = 1, 
    #          out = celltype_path )
    
  }, error = function(e) {
    # Print error message and skip to the next iteration
    message(paste("Error with celltype:", celltype, "-", e$message))
    # remove the created folder
    unlink(celltype_path, recursive = TRUE)
  })
}


