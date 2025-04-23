
# PseudoBulk_fn -----------------------------------------------------------
PseudoBulk_fn <- function(obj,vars){
  # obj: seurat object
  # vars: vector of variables to aggregate by
  # return: a list of data frames with the aggregated counts
  
  # pseudo bulk
  pseudo_obj <- AggregateExpression(obj, assays = "RNA",
                                    group.by = vars,
                                    slot = "counts",
                                    return.seurat = F)
  
  # extracting the object and making into a df
  cts <- data.frame(pseudo_obj$RNA)
  # transpose
  cts.t <- t(cts)
  # get the values where to split 
  splitRows <- gsub('_.*', '', rownames(cts.t))
  # split data.frame
  cts.split <- split.data.frame(cts.t, f = factor(splitRows))
  # fix colnames and transpose
  cts.split.modified <- lapply(cts.split, function(x){
    # Remove the first part of the colnames (cell type names)
    rownames(x)  <- gsub('^[^_]*_(.*)', '\\1', rownames(x))
    # Transpose
    t(x)
  })
  
  # # Changing the names to the original format
  # names(cts.split.modified) <- lapply(names(cts.split.modified), function(name) {
  #   name <- gsub("\\.", " ", name)           # Replace dots with spaces
  #   name <- gsub(" \\(.*\\)", "", name)      # Remove anything in parentheses
  #   name <- gsub("  ", " ", name)            # Replace double spaces with single space
  #   name <- trimws(name)                     # Trim leading and trailing whitespace
  #   return(name)
  # })
  
  return(cts.split.modified)

}

# Prepare DE table for searchlight ---------------------------------------------
searchlight_prep <- function(res){
  
  # converting to a df
  res_df <- data.frame(res)
  
  # order by padj
  res_df <- res_df[order(res_df$padj),]
  # gene(ID) column
  res_df$ID <- rownames(res_df)
  # Move 'ID' column to the first position
  res_df <- res_df[, c("ID", setdiff(names(res_df), "ID"))]
  
  # subset columns
  res_df <- select(res_df, ID, log2FoldChange, pvalue, padj)
  # rename the columns
  colnames(res_df) <- c("ID", "log2Fold", "P", "P.adj")
  
  return(res_df)
  
}

# GetDE_fn (differential expression) ------------------------------------------

GetDE_fn <- function(em_df, var, suffix, sample_sheet, reference, contrast, p_threshold, fold_threshold, out){
  # em_list: list of em df
  # var : variable(celltype/region ...)  name to be used in the output file name
  # sample_sheet: sample sheet for DESeq2
  # suffix: suffix to add to the file name
  # reference <- reference var
  # contrast <- contrast var
  # fold_threshold: default 1 
  # out: output directory where the CSVs will be saved
  
  suppressMessages({
    # CODE
    
    # subsetting only samples in ss
    counts <- em_df[, sample_sheet$sample]
    
    # filtering 
    counts <- subset(counts,apply(counts, 1, mean) >= 1)
    counts <-  as.matrix(counts)
    
    # other way of filtering
    # keep <- rowSums(counts(dds) >=10) >= 3
    # dds <- dds[keep,]
    
    # Deseq2 analysis
    dds <- DESeqDataSetFromMatrix(countData = counts,
                                  colData = sample_sheet,
                                  design = ~ sample_group)
    
    # setting  reference condition
    dds$sample_group <- relevel(dds$sample_group, ref = reference)
    
    # uncomment in case you want to check the reference and contrast
    print(paste0("Reference: ",levels(dds$sample_group)[1]))
    print(paste0("Contrast: ",levels(dds$sample_group)[2]))
    
    
    
    # run DESeq2
    dds <- DESeq(dds)
    
    
    # check dispersion shrinkage
    # plotDispEsts(dds, main="Gene Matrix", cex.lab = 1, cex.main=2)
    
    # Generate results object USING LFC
    resLFC <- lfcShrink(dds,
                        coef = 2,
                        type = "ashr", # ashr or apeglm
                        lfcThreshold = fold_threshold,
                        alpha = p_threshold)
    
    
    # converting to a df
    res_df <- data.frame(resLFC)
    
    # order by padj
    res_df <- res_df[order(res_df$padj),]
    # gene(ID) column
    res_df$ID <- rownames(res_df)
    # Move 'ID' column to the first position
    res_df <- res_df[, c("ID", setdiff(names(res_df), "ID"))]
    
    # subset columns
    res_df <- select(res_df, ID, log2FoldChange, pvalue, padj)
    # rename the columns
    colnames(res_df) <- c("ID", "log2Fold", "P", "P.adj")
    
    # write csv to out location
    write.table(res_df,
                file = file.path(out, paste0("de_", var, "_", suffix, ".tsv")),
                sep = "\t",
                row.names = FALSE,
                quote = FALSE)
    
    print(paste("DE table for ", var," | ", suffix, "saved"))
    
  })
  
}

# NOT YET NEEDED ----------------------------------------------------------

# EM table generator function ---------------------------------------------
GetEM_fn <- function(pseudo_obj, Sample){
  
  # create directory for em tables 
  create_or_overwrite_directory(file.path(output_dir, Sample, "EM"))
  create_or_overwrite_directory(file.path(output_dir, Sample, "EM_Scaled"))
  
  # Create a file to log errors
  error_log_file <- file.path(output_dir, Sample, "EM", "error_log.txt")
  cat("", file = error_log_file)  # Create or overwrite the error log file
  
  # loop to get all em tables
  for (cluster in names(pseudo_obj)) {
    tryCatch({
      
      em <- as.data.frame(pseudo_obj[[cluster]])
      
      em.scaled <- na.omit(as.data.frame(t(scale(t(em))))) # scaling the data
      em.scaled$gene <- rownames(em) # adding gene column
      
      write.csv(em.scaled, file = file.path(output_dir, Sample, "EM_Scaled", paste0("bulk_", cluster, "_em.scaled.csv")), row.names = FALSE)
      print(paste("EM scaled table for Cluster", cluster, "saved"))
      
      print(class(em))
      
      em$em_mean <- rowMeans(em) # em mean
      em$gene <- rownames(em) # adding gene column
      
      write.csv(em, file = file.path(output_dir, Sample, "EM", paste0("bulk_", cluster, "_em.csv")), row.names = FALSE)
      print(paste("EM table for Cluster", cluster, "saved"))
    }, error = function(e) {
      error_message <- paste("Error with Cluster", cluster, ":", e$message, "\n")
      cat(error_message, file = error_log_file, append = TRUE)
      print(error_message)
    })
  }
  
  final_statement <- paste("EM table generated for", Sample, "finished succesfully", sep = " ")
  
  
}

##### Master table generator function #####
GetMaster_fn <- function(pseudo_obj, Sample){
  
  # create directory for Master tables 
  create_or_overwrite_directory(file.path(output_dir, Sample, "Master"))
  
  # Create a file to log errors
  error_log_file <- file.path(output_dir, Sample, "Master", "error_log.txt")
  cat("", file = error_log_file)  # Create or overwrite the error log file
  
  # loop to merge all de and em tables 
  for (cluster in names(pseudo_obj)) {
    tryCatch({
      em <- read.csv(file.path(output_dir, Sample, "EM", paste0("bulk_", cluster, "_em.csv")))
      de <- read.csv(file.path(output_dir, Sample, "DE", paste0("bulk_", cluster, "_de.csv")))
      master <- merge(em, de, by = "gene")
      write.csv(master, file = file.path(output_dir, Sample, "Master", paste0("bulk_", cluster, "_master.csv")), row.names = FALSE)
      print(paste("Master table for Cluster", cluster, "saved"))
    }, error = function(e) {
      error_message <- paste("Error with Cluster", cluster, ":", e$message, "\n")
      cat(error_message, file = error_log_file, append = TRUE)
      print(error_message)
    })
  }
  
  final_statement <- paste("Master tables generated for", Sample, "finished succesfully", sep = " ")
  
}

##### Volcano plot generator function #####
GetVolcano_fn <- function(pseudo_obj, Sample, p_threshold, fold_threshold){
  
  # create dir for volcano plots
  create_or_overwrite_directory(file.path(output_dir, Sample, "Volcano"))
  
  # Create a file to log errors
  error_log_file <- file.path(output_dir, Sample, "Volcano", "error_log.txt")
  cat("", file = error_log_file)  # Create or overwrite the error log file
  
  # loop to get all the volcano plots 
  for (cluster in names(pseudo_obj)){
    tryCatch({
      cluster_master <- read.table(file.path(output_dir, Sample, "Master", paste0("bulk_", cluster, "_master.csv")), header = TRUE, sep = ",")
      
      volcano_pl <- volcano_plot(master = cluster_master,
                                 p_threshold = p_threshold,
                                 fold_threshold = fold_threshold,
                                 title = cluster)
      
      print(paste("Volcano Plot Cluster", cluster, "saved", sep = " "))
      
      # print hte names of the most significant genes in master
      print(head(cluster_master[order(cluster_master$p.adj),], 10)$gene)
      
      # Saving the plot
      save_plot(plot = volcano_pl,
                file_name = file.path(output_dir, Sample, "Volcano", paste0("volcano_", cluster, ".png")),
                w=3000, h=1500)
    }, error = function(e) {
      error_message <- paste("Error with Cluster", cluster, ":", e$message, "\n")
      cat(error_message, file = error_log_file, append = TRUE)
      print(error_message)
    })
  }
  
}

##### Container function #####
DiffExp_fn <- function(Sample, Ref_Cond, Contr_Cond, p_threshold, fold_threshold){
  
  # object name
  sample_labelled.rds <- paste0(Sample,"_labelled.rds")
  # load object
  obj <- readRDS(file = file.path(output_dir, Sample, sample_labelled.rds))
  
  # 1. pseudobulk
  PseudoBulk_fn(obj, Sample)
  
  #psudobulked object name
  pseudo_object.rds <- paste0("pseudo_", Sample, ".rds")
  # load the psudobulked object
  pseudo_obj <- readRDS(file = file.path(output_dir, Sample, pseudo_object.rds))
  
  # 2. differential expression
  GetDE_fn(pseudo_obj, Sample, Ref_Cond, Contr_Cond, p_threshold, fold_threshold)
  
  # 3. EM tables
  GetEM_fn(pseudo_obj, Sample)
  
  # # 4. Master tables
  GetMaster_fn(pseudo_obj, Sample)
  
  
  # 5. Volcano plots
  GetVolcano_fn(pseudo_obj, Sample, p_threshold, fold_threshold)
  
  print("******************************************")
  
}