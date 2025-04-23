
# load the broad object ---------------------------------------------------

broad <- readRDS(file.path(objects_dir,"broad.label_objects",
                           "2024_12_15_broad_pos.RDS"))


# names for output files --------------------------------------------------

pca_obj_name <- "2024_12_15_broad_pos_pca_dim_reduction.RDS"
unintegrated_obj_name <- "2024_12_15_broad_pos_unintegrated.RDS"

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

saveRDS(broad, file.path(objects_dir,"broad.label_objects", pca_obj_name))

# parameters for UMAP/clustering -----------------------------------------------
ElbowPlot(cluster.obj)
dims <- c(1:10)
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

saveRDS(broad, file.path(objects_dir,"broad.label_objects", unintegrated_obj_name ))


