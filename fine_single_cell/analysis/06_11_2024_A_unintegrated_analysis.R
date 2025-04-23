
# load the broad object ---------------------------------------------------

broad <- readRDS(file.path(objects_dir,"broad.label_objects",
                           "05_11_2024_broad_v1.0.RDS"))


# removing dimensional reductions -----------------------------------------

for (cluster in names(broad)){
  
  # indicate cluster in console
  print(cluster)
  
  cluster.obj <- broad[[cluster]]
  
  # remove reductions 
  cluster.obj@reductions <- list()
  
  # save
  broad[[cluster]] <- cluster.obj
}


# standard processing -----------------------------------------------------

for (cluster in names(broad)){
  
  # indicate cluster in console
  print(cluster)
  
  cluster.obj <- broad[[cluster]]
  
  # normalisation
  cluster.obj <- NormalizeData(cluster.obj)
  cluster.obj <- ScaleData(cluster.obj, features = rownames(cluster.obj))
  cluster.obj <- RunPCA(cluster.obj, features = rownames(cluster.obj))

  # save
  broad[[cluster]] <- cluster.obj
  
}

#  save broad object (PCA) -----------------------------------------------------

saveRDS(broad, file.path(objects_dir,"broad.label_objects",
                         "06_11_2024_broad_v1.1_pca_dim_reduction.RDS"))

# parameters for UMAP/clustering -----------------------------------------------

dims <- c(1:20)
res <- 0.2

# UMAP dimensional reduction  ---------------------------------------------

for (cluster in names(broad)){
  
  # indicate cluster in console
  print(cluster)
  
  cluster.obj <- broad[[cluster]]
  
  # UMAP
  cluster.obj <- RunUMAP(cluster.obj, dims = dims,
                         reduction = "pca",
                         reduction.name = "umap.unintegrated" )
  
  # save
  broad[[cluster]] <- cluster.obj
  
}


# clustering --------------------------------------------------------------

for (cluster in names(broad)){
  
  # indicate cluster in console
  print(cluster)
  
  cluster.obj <- broad[[cluster]]
  
  # clustering
  cluster.obj <- FindNeighbors(cluster.obj, dims = dims)
  cluster.obj <- FindClusters(cluster.obj, resolution = res, cluster.name = "unintegrated_subclusters")
  
  # save
  broad[[cluster]] <- cluster.obj
  
}


# save the broad object (clustering) -------------------------------------------

saveRDS(broad, file.path(objects_dir,"broad.label_objects",
                         "06_11_2024_broad_v1.2_unintegrated.RDS"))


