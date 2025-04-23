# load brain/regions_bulk objects
regions_bulk <- readRDS(file.path(data_dir,"regional",
                                    "regions.major_pseudo_bulk_roi.RDS"))

# creating a ss list for all cell types -------------------------------------------------
ss_list <- list()

for (region in names(regions_bulk)){
  
  # get the region
  counts <- regions_bulk[[region]]
  
  # generate sample level metadata
  ss <- tibble(sample = colnames(counts))
  
  # generating mouse_id column
  segments <- c(1) # mouse_id
  replacement <- paste0("\\", segments[1])
  ss$mouse_id <- gsub('^([^_]+)_([^_]+)+_([^_]+)$',
                      replacement, 
                      colnames(counts))
  
  # generating position column
  segments <- c(2) # position
  replacement <- paste0("\\", segments[1])
  ss$position <- gsub('^([^_]+)_([^_]+)+_([^_]+)$',
                      replacement, 
                      colnames(counts))
  
  # generating condition column
  segments <- c(3) # condition
  replacement <- paste0("\\", segments[1])
  ss$condition <- gsub('^([^_]+)_([^_]+)+_([^_]+)$',
                       replacement, 
                       colnames(counts))
  
  # generating sample group column
  ss$sample_group <- paste0(ss$position, "_", ss$condition)
  
  # print to check
  print(ss)
  
  ss_list[[region]] <- ss
  
}

# Subset ss to include only the ALDARA samples -------------------------------
for (region in names(regions_bulk)){
  
  # get the region
  ss <- ss_list[[region]]
  
  # subset to include only ALDARA samples
  ss_c <- ss[ss$condition == "A",]
  
  # print to check
  print(ss_c)
  
  # add to the ss_list
  ss_list[[region]] <- ss_c
  
}

# creating a list of em tables for all regions --------------------------------
em_list <- list()

for (region in names(regions_bulk)){
  
  cluster_df <- as.data.frame(regions_bulk[[region]])
  
  cluster_df$ID <- rownames(cluster_df)
  
  cluster_df <- cluster_df[, c("ID", setdiff(names(cluster_df), "ID"))]
  
  # add to the em_list
  em_list[[region]] <- cluster_df
  
}

# subset list of em tables to include only anterior samples ----------------------
for (region in names(regions_bulk)){
  
  # get the region
  em <- em_list[[region]]
  
  # c samples 
  a_samples <- ss_list[[region]]$sample
  
  # subset to include only anterior samples
  em_a <- em[, c("ID", a_samples)]
  
  # print to check
  print(em_a)
  
  # add to the em_list
  em_list[[region]] <- em_a
  
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
    base_path <- file.path(data_dir,"regional","differential_expression", "position", "aldara")
    
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
    # subset ss to include only 
    # save the ss
    write.table(ss,
                file = file.path(region_path, paste0("ss_", region, ".tsv")),
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
    
    ## DESeq
    
    ## Create DESeq object
    dds <- DESeqDataSetFromMatrix(countData = em[,ss$sample], 
                                  colData = ss,
                                  design = ~ mouse_id + sample_group)
    
    # Filter lowly expressed genes
    dds <- dds[rowMeans(counts(dds)) >= 1, ] # Keep genes with mean raw counts >= 1 
    
    # setting  reference condition
    dds$sample_group <- relevel(dds$sample_group, ref = "ant_A")
    
    ## Run DESEQ2
    dds <- DESeq(dds)
    
    #Checking coeficients
    resultsNames(dds)
    
    ## results
    res <- lfcShrink(dds,
                     coef = 2,
                     type = "ashr") # ashr or apeglm
    
    res_df <- searchlight_prep(res)
    
    # save the results
    write.table(res_df,
                file = file.path(region_path, paste0("de_", region, ".tsv")),
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
    
    
  }, error = function(e) {
    # Print error message and skip to the next iteration
    message(paste("Error with region:", region, "-", e$message))
    # remove the created folder
    unlink(region_path, recursive = TRUE)
  })
}


