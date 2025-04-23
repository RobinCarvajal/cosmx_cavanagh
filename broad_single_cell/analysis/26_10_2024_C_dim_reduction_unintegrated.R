# load qc object
brain <- readRDS(file.path(main_dir, "New_Analysis","objects",
                              "25_10_2024_brain_v1.1_QC.RDS"))

# normalisation and PCA reduction --------------------------------------------

n_dims <- 10 # Megan used 10 dims
res <- 0.7 # Megan used 0.7 resolution
n_var_features <- 1000 # I used 1000 varaible features (all)

brain <- NormalizeData(brain)
brain <- FindVariableFeatures(brain, nfeatures = n_var_features)
brain <- ScaleData(brain)
brain <- RunPCA(brain)

# elbow plot
ElbowPlot(brain, ndims = 50) 


# Find Neighbours/Clusters ------------------------------------------------

brain <- FindNeighbors(brain, dims = 1:n_dims, reduction = "pca")
brain <- FindClusters(brain, resolution = res,
                         cluster.name = "unintegrated_clusters")

# UMAP unintegrated -------------------------------------------------------

brain <- RunUMAP(brain, dims = 1:n_dims, reduction = "pca", reduction.name = "umap.unintegrated")

# Dimplot
DimPlot(brain, reduction = "umap.unintegrated", group.by = "unintegrated_clusters")
DimPlot(brain, reduction = "umap.unintegrated", group.by = "condition")

# save the object 
saveRDS(brain, file = file.path(main_dir, "New_Analysis","objects",
                                   "26_10_2024_brain_v1.2_unintegrated.RDS"))
