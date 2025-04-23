
#read pseudo object 
pseudo_brain <- readRDS(file.path(data_dir,"single_cell","differential_expression",
                                  "30_10_2024_pseudo_bulk.RDS"))


# Code  --------------------------------------------------------------------

for (celltype in names(pseudo_brain)){
  
  # out path 
  celltype_path <- file.path(data_dir,"single_cell","differential_expression", "cell_types",celltype)
  # make dir
  if (!dir.exists(celltype_path)) {dir.create(celltype_path, recursive = TRUE)}
  
  cluster_df <- as.data.frame(pseudo_brain[[celltype]])
  
  cluster_df$ID <- rownames(cluster_df)
  
  cluster_df <- cluster_df[, c("ID", setdiff(names(cluster_df), "ID"))]
  
  # save 
  write.table(cluster_df, 
              file = file.path(celltype_path, paste0("em_", celltype, ".tsv")), 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)
  
}