# load fine pseudobulk object 
pseudo_brain <- readRDS(file.path(data_dir,"de_utils",
                                  "fine.labels_pseudo_bulk.RDS"))

all_samples <- colnames(pseudo_brain[["Astrocytes.0"]])

test <- as.data.frame(pseudo_brain[["Astrocytes.0"]])
ncol(test)


failed <- c()
for (celltype in names(pseudo_brain)){
  test <- as.data.frame(pseudo_brain[[celltype]])
  if (ncol(test) < 16){
    print(celltype)
    print(ncol(test))
    failed <- c(failed, celltype)
  }
}
failed

failed_list <- list()
for (celltype in failed){
  test <- as.data.frame(pseudo_brain[[celltype]])
  # difference between all samples and the samples in each froup
  missing_samples <- setdiff(all_samples, colnames(test))
  test_samples <- list(celltype = celltype, 
                       #samples = colnames(test),
                       n_missing = 16 - ncol(test),
                       missing_samples = missing_samples)
  failed_list[[celltype]] <- test_samples
}
failed_list

