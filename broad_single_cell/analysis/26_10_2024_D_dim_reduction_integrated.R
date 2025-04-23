# load un-integrated object
brain <- readRDS(file.path(main_dir, "New_Analysis","objects",
                           "26_10_2024_brain_v1.2_unintegrated.RDS"))


# parameters --------------------------------------------------------------

# harmony
harmony_arguments <- c("condition", "sample_name") # Megan controlled for sample only
harmony_theta <- c(1,1)
harmony_red_name <- "harmony"

# clustering
clusters_name <- "integrated_clusters"
n_dims <- 10 
res <- 0.7 

# UMAP
umap_red_name <- "umap.integrated"

# Harmony PCA embeddings batch correction ---------------------------------
brain <- RunHarmony(brain,
                    group.by.vars = harmony_arguments, 
                    theta = harmony_theta,
                    reduction = "pca", 
                    assay.use = "RNA", 
                    reduction.save = harmony_red_name,
                    verbose=FALSE)

# Find Neighbours/Clusters ------------------------------------------------

brain <- FindNeighbors(brain, dims = 1:n_dims, reduction = harmony_red_name)
brain <- FindClusters(brain, resolution = res,
                         cluster.name = clusters_name)

# UMAP integrated -------------------------------------------------------

brain <- RunUMAP(brain, dims = 1:n_dims, 
                 reduction = harmony_red_name, 
                 reduction.name = umap_red_name)


# save RDS
saveRDS(brain, file.path(main_dir, "New_Analysis","objects",
                         "26_10_2024_brain_v1.3_integrated.RDS"))
