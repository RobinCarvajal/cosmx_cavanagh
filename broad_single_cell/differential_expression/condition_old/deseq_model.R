# load brain/regions_bulk objects
regions_bulk <- readRDS(file.path(data_dir,"regional",
                                  "2025_01_23_region-bulk.RDS"))

# creating a list of em tables for all regions --------------------------------
em_list <- list()

for (region in names(regions_bulk)){
  
  cluster_df <- as.data.frame(regions_bulk[[region]])
  
  cluster_df$ID <- rownames(cluster_df)
  
  cluster_df <- cluster_df[, c("ID", setdiff(names(cluster_df), "ID"))]
  
  # add to the em_list
  em_list[[region]] <- cluster_df
  
}

# creating a ss list for all cell types -------------------------------------------------
ss_list <- list()

for (region in names(regions_bulk)){
  
  # get the region
  counts <- regions_bulk[[region]]
  
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
  
  ss_list[[region]] <- ss
  
}

# ALL comparisons  -------------------------------------------------------

# NOTE: I do all comparisons in a loop to avoid repeating code

# in a loop
for (region in names(regions_bulk)) {
  
  # current region
  print(region)
  
  # Wrap in tryCatch to handle errors
  tryCatch({
    # base path
    base_path <- file.path(data_dir,"regional","differential_expression", "condition")
    
    # out path 
    region_path <- file.path(base_path, region)
    # make dir
    if (!dir.exists(region_path)) {dir.create(region_path, recursive = TRUE)}
    
    # em 
    em <- em_list[[region]]
    # save the em 
    write.table(em, 
                file = file.path(region_path, paste0("em_", region, ".tsv")), 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
    # ss table
    ss <- ss_list[[region]]
    # save the ss
    write.table(ss,
                file = file.path(region_path, paste0("ss_", region, ".tsv")),
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
    
    # Control vs Aldara
    GetDE_fn(em_df = em,
             var = region,
             suffix = "C_vs_A", 
             sample_sheet = ss,
             reference = "C", 
             contrast = "A", 
             p_threshold = 0.5, 
             fold_threshold = 1, 
             out = region_path )
    
  }, error = function(e) {
    # Print error message and skip to the next iteration
    message(paste("Error with region:", region, "-", e$message))
    # remove the created folder
    unlink(region_path, recursive = TRUE)
  })
}
