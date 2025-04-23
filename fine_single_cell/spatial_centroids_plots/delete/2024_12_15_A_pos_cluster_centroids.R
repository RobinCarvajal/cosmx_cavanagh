# read broad object
broad <- readRDS(file.path(objects_dir,"broad.label_objects",
                           "2024_12_15_broad_pos_integrated.RDS"))

# create a plots directory 
plots_dir <- file.path(data_dir, "fine_single_cell", "pos", "spatial_centroids_plots")

# Create the directory if they do not already exist
if (!dir.exists(plots_dir)) {dir.create(plots_dir, recursive = TRUE)}


# centroids plots -------------------------------------------------------------

for (cluster in names(broad)){
  
  cluster.obj <- broad[[cluster]]
  
  fov_list <- names(cluster.obj@images)
  
  for (fov in fov_list){
    
    
    # plot the centroids
    fov_plot <- ImageDimPlot(cluster.obj, fov = fov, 
                             split.by = "integrated_subclusters",
                             flip_xy = FALSE)+
      ggtitle(paste0(cluster," : ", fov)) +
      theme(plot.title = element_text(size = 40, hjust = 0.5))
    
    # create a folder plot per cluster
    cluster_plots_dir <- file.path(plots_dir, cluster)
    if (!dir.exists(cluster_plots_dir)) {dir.create(cluster_plots_dir, recursive = TRUE)}
    
    
    # save the plots
    ggsave(plot = fov_plot,
           file.path(cluster_plots_dir, paste0(fov, "_centroids.png")),
           width = 20, height = 20, dpi = 600)
    
  }
  
}