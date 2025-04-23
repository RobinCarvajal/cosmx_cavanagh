
# load data
brain <- readRDS(file.path(data_dir, "objects","brain_objects",
                           "2025_01_22_brain_v1.7_region.labels_added.RDS"))

results_dir <- file.path("F:/projects_work/cosmx_cavanagh/cosmx_results/requests/new_umap")
if (!dir.exists(results_dir)) {dir.create(results_dir, recursive = TRUE)}

# test2 - same a test 1 bus using uwot.sgd=TRUE should be faster
Sys.time()
brain <- RunUMAP(brain, dims = 1:10, 
                 reduction = "harmony", 
                 reduction.name = "umap.integrated",
                 uwot.sgd = TRUE, 
                 n.neighbors = 30, 
                 min.dist = 0.05 , # (distance within clusters)
                 spread = 0.4) # (distance between clusters)
DimPlot(brain, reduction = "umap.integrated", group.by = "integrated_clusters", label.box = T, label = T)
Sys.time()


png(file.path(results_dir, "2025_01_22_brain_v1.7_region.labels_added.png"),
    width = 10, height = 10, res = 600, units="in")
DimPlot(brain, reduction = "umap.integrated", group.by = "integrated_clusters", label.box = T, label = T)
dev.off()

# save the object
svg(file.path(results_dir, "2025_01_22_brain_v1.7_region.labels_added.svg"), width = 10, height = 10)
DimPlot(brain, reduction = "umap.integrated", group.by = "integrated_clusters", label.box = T, label = T)
dev.off()

# save object 
saveRDS(brain, file.path(data_dir, "objects","brain_objects",
                           "2025_01_22_brain_v1.7_region.labels_added_new_umap.RDS"))


