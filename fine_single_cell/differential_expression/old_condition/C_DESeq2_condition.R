# load brain/pseudo_brain objects
pseudo_brain <- readRDS(file.path(objects_dir,"pseudo_objects",
                                  "fine.labels_pseudo_bulk.RDS"))

# creating a ss list for all cell types -------------------------------------------------
ss_list <- list()

for (celltype in names(pseudo_brain)){
  
  # get the celltype
  counts <- pseudo_brain[[celltype]]
  
  # generate sample level metadata
  ss <- tibble(sample = colnames(counts))
  
  segments <- c(4) # condition
  
  replacement <- paste0("\\", segments[1])
  
  # modified so that Searchligth2 will accept it
  ss$sample_group <- gsub('^([^_]+)_([^_]+)+_([^_]+)_([^_]+)$',
                          replacement, 
                          colnames(counts))
  
  # print to check
  print(ss)
  
  ss_list[[celltype]] <- ss
  
}

# ALL comparisons  -------------------------------------------------------

# NOTE: I do all comparisons in a loop to avoid repeating code

# in a loop
for (celltype in names(pseudo_brain)) {
  
  # current celltype
  print(celltype)
  
  # Wrap in tryCatch to handle errors
  tryCatch({
    # base path
    base_path <- file.path(data_dir,"fine_single_cell","differential_expression", "condition", "cell_types")
    
    # out path 
    celltype_path <- file.path(base_path, celltype)
    # make dir
    if (!dir.exists(celltype_path)) {dir.create(celltype_path, recursive = TRUE)}
    
    # em path 
    em <- pseudo_brain[[celltype]]
    
    ss <- ss_list[[celltype]]
    
    # Control vs Aldara
    GetDE_fn(em_df = em,
             celltype = celltype,
             suffix = "C_vs_A", 
             sample_sheet = ss,
             reference = "C", 
             contrast = "A", 
             p_threshold = 0.5, 
             fold_threshold = 1, 
             out = celltype_path )
  }, error = function(e) {
    # Print error message and skip to the next iteration
    message(paste("Error with celltype:", celltype, "-", e$message))
    # remove the created folder
    unlink(celltype_path, recursive = TRUE)
  })
}
