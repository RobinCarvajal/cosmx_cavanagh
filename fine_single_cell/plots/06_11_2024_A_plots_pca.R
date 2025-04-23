broad <- readRDS(file.path(objects_dir,"broad.label_objects",
                           "06_11_2024_broad_v1.1_pca_dim_reduction.RDS"))

# create a plots directory 
plots_dir <- file.path(data_dir, "broad_single_cell", "plots")

# Create the directory if they do not already exist
if (!dir.exists(plots_dir)) {dir.create(plots_dir, recursive = TRUE)}


# pca plots ---------------------------------------------------------------

for (cluster in names(broad)){
  
  cluster.obj <- broad[[cluster]]
  
  # elbow plot
  plot1 <- ElbowPlot(cluster.obj, ndims = 50)
  # heatmap dims
  plot2 <- DimHeatmap(cluster.obj, dims = 1:20, cells = 500, balanced = TRUE, fast = FALSE)
  
  # create a folder plot per cluster
  cluster_plots_dir <- file.path(plots_dir, cluster)
  if (!dir.exists(cluster_plots_dir)) {dir.create(cluster_plots_dir, recursive = TRUE)}
  
  # save the plots
  ggsave(plot = plot1,
         file.path(cluster_plots_dir, "elbow_plot.png"),
         width = 10, height = 5, dpi = 300, bg = "white")
  
  ggsave(plot = plot2,
         file.path(cluster_plots_dir, "heatmap_dims.png"),
         width = 7, height = 10, dpi = 300)
  
}
