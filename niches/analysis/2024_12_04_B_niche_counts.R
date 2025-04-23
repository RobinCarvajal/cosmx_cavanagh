library(gridExtra)
# Load niches list --------------------------------------------------------
ant_niches.list <- readRDS(file.path(data_dir,"niches","2024_12_04_niches.list_ant_k5.RDS"))

pos_niches.list <- readRDS(file.path(data_dir,"niches","2024_12_04_niches.list_pos_k_8_12.RDS"))


# ANTERIOR niches --------------------------------------------------------

for (sample in names(ant_niches.list)){
  
  # sample object
  sample_obj <- ant_niches.list[[sample]]@meta.data
  
  # niches k 
  niches_k <- paste0("niches_k",5)
  
  # Calculating Cell counts
  cell_count <- table(sample_obj[["broad.labels"]], sample_obj[[niches_k]])
  # Calculate row totals and add as a new column
  cell_count <- cbind(cell_count, Column_Totals = rowSums(cell_count))
  # Calculate column totals and add as a new row
  cell_count <- rbind(cell_count, Row_Totals = colSums(cell_count))
  
  # Define file path for the PDF
  pdf_file1 <- file.path(data_dir, "niches", sample, paste0(sample, "_cell_counts.pdf"))
  # Create the PDF
  pdf(pdf_file1, width = 8.5, height = 11) # Adjust dimensions as needed
  # Convert the data frame into a table format
  grid.table(cell_count)
  # Close the PDF device
  dev.off()
  
  # add celltype column 
  cell_count <- cbind(cell_count, celltype = rownames(cell_count))
  # set as 1st column
  cell_count <- cell_count[,c(ncol(cell_count),1:ncol(cell_count)-1)]
  
  # save table as tsv
  write.table(cell_count, 
              file.path(data_dir,"niches",sample, paste0(sample,"_cell_count.tsv")), 
              sep = "\t",
              quote=F,
              row.names=F)
  
  ## Proportions 
  cell_proportions <- table(sample_obj[["broad.labels"]], sample_obj[[niches_k]])
  cell_proportions <- prop.table(cell_proportions, margin = 2)
  cell_proportions <- round(cell_proportions*100,2)
  
  
  # Define file path for the PDF
  pdf_file2 <- file.path(data_dir, "niches", sample, paste0(sample, "_cell_proportions.pdf"))
  # Create the PDF
  pdf(pdf_file2, width = 8.5, height = 11) # Adjust dimensions as needed
  # Convert the data frame into a table format
  grid.table(cell_proportions)
  # Close the PDF device
  dev.off()
  
}


# note also do proportions  and other counts...



# POSTERIOR niches --------------------------------------------------------


for (sample in names(pos_niches.list)){
  
  # sample object metadata
  sample_obj <- pos_niches.list[[sample]]@meta.data
  
  # niches k 
  niches_k <- paste0("niches_k",7)
  
  # Calculating Cell counts
  cell_count <- table(sample_obj[["broad.labels"]], sample_obj[[niches_k]])
  # Calculate row totals and add as a new column
  cell_count <- cbind(cell_count, Column_Totals = rowSums(cell_count))
  # Calculate column totals and add as a new row
  cell_count <- rbind(cell_count, Row_Totals = colSums(cell_count))
  
  # Define file path for the PDF
  pdf_file1 <- file.path(data_dir, "niches", sample, paste0(sample, "_cell_counts.pdf"))
  # Create the PDF
  pdf(pdf_file1, width = 8.5, height = 11) # Adjust dimensions as needed
  # Convert the data frame into a table format
  grid.table(cell_count)
  # Close the PDF device
  dev.off()
  
  # add celltype column 
  cell_count <- cbind(cell_count, celltype = rownames(cell_count))
  # set as 1st column
  cell_count <- cell_count[,c(ncol(cell_count),1:ncol(cell_count)-1)]
  
  # save table as tsv
  write.table(cell_count, 
              file.path(data_dir,"niches",sample, paste0(sample,"_cell_count.tsv")), 
              sep = "\t",
              quote=F,
              row.names=F)
  
  ## Proportions 
  cell_proportions <- table(sample_obj$broad.labels, sample_obj[[niches_k]])
  cell_proportions <- prop.table(cell_proportions, margin = 2)
  cell_proportions <- round(cell_proportions*100,2)
  
  
  # Define file path for the PDF
  pdf_file2 <- file.path(data_dir, "niches", sample, paste0(sample, "_cell_proportions.pdf"))
  # Create the PDF
  pdf(pdf_file2, width = 8.5, height = 11) # Adjust dimensions as needed
  # Convert the data frame into a table format
  grid.table(cell_proportions)
  # Close the PDF device
  dev.off()
  
}