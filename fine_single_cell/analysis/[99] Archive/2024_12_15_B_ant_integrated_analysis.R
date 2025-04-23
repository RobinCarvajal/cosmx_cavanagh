
# load the broad object ---------------------------------------------------

broad <- readRDS(file.path(objects_dir,"broad.label_objects",
                           "2024_12_15_broad_ant_unintegrated.RDS"))

# names for output files --------------------------------------------------

integrated_obj_name <- "2024_12_15_broad_ant_integrated.RDS"

# integration parameters --------------------------------------------------

# harmony
harmony_arguments <- c("condition", "sample_name") # Megan controlled for sample only
harmony_theta <- c(2,2)
harmony_red_name <- "harmony"

# clustering
clusters_name <- "integrated_subclusters"
dims <- c(1:10)
res <- 0.2

# UMAP
umap_red_name <- "umap.integrated"


# Harmony PCA correction -----------------------------------------------------

for (cluster in names(broad)){
  
  # indicate cluster in console
  print(cluster)
  
  cluster.obj <- broad[[cluster]]
  
  # harmony pca correction
  cluster.obj <- RunHarmony(cluster.obj,
                      group.by.vars = harmony_arguments, 
                      theta = harmony_theta,
                      reduction = "pca", 
                      assay.use = "RNA", 
                      reduction.save = harmony_red_name,
                      verbose=FALSE)
  
  
  # save
  broad[[cluster]] <- cluster.obj
  
}

# UMAP ------------------------------------------------------------------

for (cluster in names(broad)){
  
  # indicate cluster in console
  print(cluster)
  
  cluster.obj <- broad[[cluster]]
  
  # UMAP
  cluster.obj <- RunUMAP(cluster.obj, dims = dims, 
                         reduction = harmony_red_name,
                         reduction.name = umap_red_name)
  
  # save
  broad[[cluster]] <- cluster.obj
  
}


# clustering ------------------------------------------------------------

for (cluster in names(broad)){
  
  # indicate cluster in console
  print(cluster)
  
  cluster.obj <- broad[[cluster]]
  
  # clustering
  cluster.obj <- FindNeighbors(cluster.obj, dims = dims)
  cluster.obj <- FindClusters(cluster.obj, resolution = res, cluster.name = clusters_name)
  
  # save
  broad[[cluster]] <- cluster.obj
  
}



# save the obj ------------------------------------------------------------

saveRDS(broad, file.path(objects_dir,"broad.label_objects", integrated_obj_name))















