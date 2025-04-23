# run in WSL
source("/mnt/f/CosMx_Analysis/Set_up_WSL.R")

library(future)
options(future.globals.maxSize = 8 * 1024^3)  # Set to 8 GiB

# load qc object
brain_qc <- readRDS(file.path(main_dir, "New_Analysis","objects","25_10_2024_v1.0_brain_QC.RDS"))

# normalisation and PCA reduction --------------------------------------------

n_dims <- 10 # Megan used 10 dims
res <- 0.7 # Megan used 0.7 resolution
n_var_features <- 1000 # I used 1000 varaible features (all)

brain_qc <- NormalizeData(brain_qc)
brain_qc <- FindVariableFeatures(brain_qc, nfeatures = n_var_features)
brain_qc <- ScaleData(brain_qc)
brain_qc <- RunPCA(brain_qc, npcs = 50)

# elbow plot
ElbowPlot(brain_qc, ndims = 50) 


# Find Neighbours/Clusters ------------------------------------------------

brain_qc <- FindNeighbors(brain_qc, dims = 1:n_dims, reduction = "pca", k.param = 30,
                          annoy.metric ="cosine")
brain_qc <- FindClusters(brain_qc, resolution = res,
                         cluster.name = "unintegrated_clusters")

# UMAP unintegrated -------------------------------------------------------

brain_qc <- RunUMAP(brain_qc, dims = 1:n_dims, reduction = "pca", 
                    reduction.name = "umap.unintegrated",
                    n.neighbors = 30, min.dist = 0.01, spread = 5)

# Dimplot
umap_plot <- DimPlot(brain_qc, reduction = "umap.unintegrated", group.by = "unintegrated_clusters")
DimPlot(brain_qc, reduction = "umap.unintegrated", group.by = "condition")

# save umap plot
ggsave(file.path(main_dir, "New_Analysis","plots","26_10_2024_umap_unintegrated_nanostring_settings.png"), 
       plot = umap_plot, w = 10, h = 10, dpi = 300)

# save the object 
saveRDS(brain_qc, file = file.path(main_dir, "New_Analysis","objects","26_10_2024_v1.0_brain_unintegrated_nanostring_settings.RDS"))
