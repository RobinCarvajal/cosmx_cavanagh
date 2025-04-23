
# Clustering Function -----------------------------------------------------

clustering_fn = function(obj, dims=1:10, reduction="pca", resolution=0.5, cluster.name="seurat_clusters"){
  
  # KNN graph 
  obj <- FindNeighbors(obj, reduction = "pca", dims=dims, verbose=FALSE) 
  
  # modularity optimisation: group cells together 
  obj <- FindClusters(obj, resolution=resolution, cluster.name = cluster.name, verbose=FALSE)
  
  #printing clusters and number of cells 
  print(paste0("No. clusters: ", length(unique(Idents(obj)))))
  print(table(Idents(obj)))
  
  return(obj)
  
}

