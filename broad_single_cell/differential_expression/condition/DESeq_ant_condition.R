# load brain/celltypes_bulk objects

# Immune+Astrocytes
celltypes_bulk <- readRDS(file.path(data_dir,"broad_single_cell",
                                  "broad.labels_pseudo_bulk_astro+immune.RDS"))
# Rest of cell types 
celltypes_bulk <- readRDS(file.path(data_dir,"broad_single_cell",
                                    "broad.labels_pseudo_bulk_rest.RDS"))

# creating a ss list for all cell types -------------------------------------------------
ss_list <- list()

for (celltype in names(celltypes_bulk)){
  
  # get the celltype
  counts <- celltypes_bulk[[celltype]]
  
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
  
  ss_list[[celltype]] <- ss
  
}

# Subset ss to include only the anterior samples -------------------------------
for (celltype in names(celltypes_bulk)){
  
  # get the celltype
  ss <- ss_list[[celltype]]
  
  # subset to include only anterior samples
  ss_ant <- ss[ss$position == "ant",]
  
  # print to check
  print(ss_ant)
  
  # add to the ss_list
  ss_list[[celltype]] <- ss_ant
  
}

# creating a list of em tables for all celltypes --------------------------------
em_list <- list()

for (celltype in names(celltypes_bulk)){
  
  cluster_df <- as.data.frame(celltypes_bulk[[celltype]])
  
  cluster_df$ID <- rownames(cluster_df)
  
  cluster_df <- cluster_df[, c("ID", setdiff(names(cluster_df), "ID"))]
  
  # add to the em_list
  em_list[[celltype]] <- cluster_df
  
}

# subset list of em tables to include only anterior samples ----------------------
for (celltype in names(celltypes_bulk)){
  
  # get the celltype
  em <- em_list[[celltype]]
  
  # ant samples 
  ant_samples <- ss_list[[celltype]]$sample
  
  # subset to include only anterior samples
  em_ant <- em[, c("ID", ant_samples)]
  
  # print to check
  print(em_ant)
  
  # add to the em_list
  em_list[[celltype]] <- em_ant
  
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
    base_path <- file.path(data_dir,"broad_single_cell","differential_expression", "condition", "anterior")
    
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
    # subset ss to include only 
    # save the ss
    write.table(ss,
                file = file.path(celltype_path, paste0("ss_", celltype, ".tsv")),
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")
    
    ## DESeq
    
    ## Create DESeq object
    dds <- DESeqDataSetFromMatrix(countData = em[,ss$sample], 
                                  colData = ss,
                                  design = ~ sample_group)
    
    # Filter lowly expressed genes
    dds <- dds[rowMeans(counts(dds)) >= 1, ] # Keep genes with mean raw counts >= 1 
    
    # setting  reference condition
    dds$sample_group <- relevel(dds$sample_group, ref = "ant_C")
    
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
                file = file.path(celltype_path, paste0("de_", celltype, ".tsv")),
                row.names = FALSE,
                quote = FALSE,
                sep = "\t")

    
  }, error = function(e) {
    # Print error message and skip to the next iteration
    message(paste("Error with celltype:", celltype, "-", e$message))
    # remove the created folder
    unlink(celltype_path, recursive = TRUE)
  })
}


