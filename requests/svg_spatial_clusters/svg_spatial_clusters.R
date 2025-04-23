# read the object
brain <- readRDS(file.path(data_dir, "objects","brain_objects",
                           "2025_01_22_brain_v1.7_region.labels_added.RDS"))

# create a spatial plot of each cell cluster in each image 
# using the spatial centroids

# dir to store the plots
svg_segmentation_plots_dir <- file.path("F:/projects_work/cosmx_cavanagh/cosmx_results/requests/svg_spatial_clusters")

# Create the directory if they do not already exist
if (!dir.exists(svg_segmentation_plots_dir)) {dir.create(svg_segmentation_plots_dir, recursive = TRUE)}
# all fov names 
fov_list <- c("s1055ant_bs1","s1055ant_bs4","s1055pos_bs3","s1055pos_bs1")

for (fov in fov_list){
  
  # plot the centroids
  fov_plot <- ImageDimPlot(brain, fov = fov, 
                           split.by = "integrated_clusters",
                           boundaries="segmentation",
                           flip_xy = FALSE,
                           border.size = 0.05) +
    NoLegend()
  
  # save the plot
  ggsave(plot = fov_plot, 
         filename = file.path(svg_segmentation_plots_dir, paste0("28_10_2024",fov,"_centroids.png")),
         w=20, h=20, dpi=600)
  
  # save as svg
  svg(file.path(svg_segmentation_plots_dir, paste0(fov,".svg")),
      width = 20, height = 20)
  print(fov_plot)
  dev.off()
  
}

