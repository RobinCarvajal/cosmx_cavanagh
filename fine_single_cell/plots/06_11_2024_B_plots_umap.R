broad <- readRDS(file.path(objects_dir,"broad.label_objects",
                           "06_11_2024_broad_v1.3_integrated.RDS"))

# create a plots directory 
plots_dir <- file.path(data_dir, "broad_single_cell", "plots")

# Create the directory if they do not already exist
if (!dir.exists(plots_dir)) {dir.create(plots_dir, recursive = TRUE)}


# UMAP unintegrated plots ---------------------------------------------------------------

for (cluster in names(broad)){
  
  # indicate cluster in console
  print(cluster)
  
  cluster.obj <- broad[[cluster]]
  
  # umap plot
  plot1 <- DimPlot(cluster.obj, reduction = "umap.unintegrated", label = TRUE, label.box = T,  
                   pt.size = 0.5, repel = TRUE, group.by = "unintegrated_subclusters")
  
  # create a folder plot per cluster
  cluster_plots_dir <- file.path(plots_dir, cluster)
  if (!dir.exists(cluster_plots_dir)) {dir.create(cluster_plots_dir, recursive = TRUE)}
  
  # save the plots
  ggsave(plot = plot1,
         file.path(cluster_plots_dir, paste0(cluster, "_unintegrated_umap.png")),
         width = 10, height = 10, dpi = 300)
  
}

# UMAP integrated --------------------------------------------------------------

for (cluster in names(broad)){
  
  # indicate cluster in console
  print(cluster)
  
  cluster.obj <- broad[[cluster]]
  
  # umap plot
  plot1 <- DimPlot(cluster.obj, reduction = "umap.integrated", label = TRUE, label.box = T,  
                   pt.size = 0.5, repel = TRUE, group.by = "integrated_subclusters")
  
  # create a folder plot per cluster
  cluster_plots_dir <- file.path(plots_dir, cluster)
  if (!dir.exists(cluster_plots_dir)) {dir.create(cluster_plots_dir, recursive = TRUE)}
  
  # save the plots
  ggsave(plot = plot1,
         file.path(cluster_plots_dir, paste0(cluster,"_integrated_umap.png")),
         width = 10, height = 10, dpi = 300)
  
}