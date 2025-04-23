
# read the broad object ---------------------------------------------------

broad <- readRDS(file.path(objects_dir,"broad.label_objects",
                           "06_11_2024_broad_v1.3_integrated.RDS"))


# read broad.markers object -----------------------------------------------

broad.markers <- readRDS(file.path(data_dir,"broad_single_cell","cluster_markers",
                                   "06_11_2024_all_broad_markers.RDS"))

# create a plots directory 
plots_dir <- file.path(data_dir, "broad_single_cell", "marker_heatmaps_plots")

# Create the directory if they do not already exist
if (!dir.exists(plots_dir)) {dir.create(plots_dir, recursive = TRUE)}

# Cluster Markers Heatmaps -----------------------------------------------------

for (cluster in names(broad)){
  
  cluster.obj <- broad[[cluster]]
  
  # downsample the object (faster plotting)
  cluster.obj <- subset(x = cluster.obj, downsample = 500)
  
  # get the markers for the cluster
  cluster.markers <-  broad.markers[[cluster]]
  
  print(cluster)
  
  # parsing the top 15 markers per cluster
  cluster.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.58) %>% # 1.5 times fold change
    slice_head(n = 15) %>%
    ungroup() -> top15
  
  # do the heatmap
  markers_heatmap <- DoHeatmap(cluster.obj, features = top15$gene) + 
    NoLegend() +
    ggtitle(cluster)+
    #center the title
    theme(
      plot.title = element_text(hjust = 0.5)
      
    )
  
  # save the heatmap
  ggsave(plot = markers_heatmap,
         filename = file.path(plots_dir, paste0(cluster,"_heatmap.png")),
         w=7, h=12, dpi=300)
  

}
