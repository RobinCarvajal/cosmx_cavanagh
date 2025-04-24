# read the object
brain <- readRDS(file.path(data_dir, "objects","brain_objects",
                           "2025_01_22_brain_v1.7_region.labels_added.RDS"))


# create a spatial plot of each cell cluster in each image 
# using the spatial centroids

# dir to store the plots
centroids_plots_dir <- file.path(main_dir, "New_Analysis", "spatial_centroids_plots")

# Create the directory if they do not already exist
if (!dir.exists(centroids_plots_dir)) {dir.create(centroids_plots_dir, recursive = TRUE)}

# all fov names 
fov_list <- unique(brain$sample_name)

for (fov in fov_list){
  
  # plot the centroids
  fov_plot <- ImageDimPlot(brain, fov = fov, 
                           split.by = "integrated_clusters",
                           flip_xy = FALSE)
  
  # save the plot
  ggsave(plot = fov_plot, 
         filename = file.path(centroids_plots_dir, paste0("28_10_2024",fov,"_centroids.png")),
         w=20, h=20, dpi=600)
}

