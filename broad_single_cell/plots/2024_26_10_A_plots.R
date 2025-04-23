
# load the integrated object
brain <- readRDS(file.path(data_dir,"objects","working_objects",
                           "26_10_2024_brain_v1.3_integrated.RDS"))

# load markers
brain.markers <- readRDS(file.path(data_dir,"single_cell","cluster_markers",
                                  "26_10_2024_all_markers.RDS"))


# create plots dir --------------------------------------------------------

# create a dir to store the plots
plots_dir <- file.path(data_dir,"single_cell", "plots")
# create a dir to store the cluster heatmaps
heatmaps_dir <- file.path(data_dir,"single_cell","cluster_heatmaps")

# Create the directory if they do not already exist
if (!dir.exists(plots_dir)) {dir.create(plots_dir, recursive = TRUE)}
if (!dir.exists(heatmaps_dir)) {dir.create(heatmaps_dir, recursive = TRUE)}


# Un-integrated UMAP -------------------------------------------------------

unintegrated_UMAP <- DimPlot(brain, 
                             reduction = "umap.unintegrated",
                             group.by = "unintegrated_clusters",
                             label.box = TRUE,
                             label = TRUE)

unintegrated_UMAP


# Integrated UMAP ---------------------------------------------------------


integrated_UMAP <- DimPlot(brain, 
                           reduction = "umap.integrated",
                           group.by = "integrated_clusters",
                           label.box = TRUE,
                           label = TRUE)

integrated_UMAP


# integrated UMAP (conditions) --------------------------------------------

integrate_UMAP_condition <- DimPlot(brain, 
                                    reduction = "umap.integrated",
                                    group.by = "condition",
                                    label.box = TRUE,
                                    label = TRUE)


# Cluster Markers Heatmaps ---------------------------------------------------------

brain.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.58) %>% # 1.5 times fold change
  slice_head(n = 30) %>%
  ungroup() -> top30

# need to downsample brain obj (faster)
brain_downsampled <- subset(x = brain, downsample = 100)

# per cluster in top20, do a heatmap
for (cluster_number in unique(top30$cluster)){
  
  # get the top 10 genes
  cluster_markers <- top30 %>% 
    dplyr::filter(cluster == cluster_number) 
  
  # do the heatmap
  markers_heatmap <- DoHeatmap(brain_downsampled, features = cluster_markers$gene) + 
    NoLegend() +
    ggtitle(paste0("Cluster ", cluster_number))
  
  # save the heatmap
  ggsave(plot = markers_heatmap,
         filename = file.path(heatmaps_dir, paste0("26_10_2024_cluster_",cluster_number,"_heatmap.png")),
         w=7, h=4, dpi=300)
  
}

# save all plots ----------------------------------------------------------

# unintegrated heatmap
ggsave(plot = unintegrated_UMAP, 
       filename = file.path(plots_dir,
                            "26_10_2024_unintegrated_UMAP.png"),
       w=10, h=10, dpi=300)

# integrated heatmap
ggsave(plot = integrated_UMAP, 
       filename = file.path(plots_dir,
                            "26_10_2024_integrated_UMAP.png"),
       w=10, h=10, dpi=300)

# integrated umap condition 
ggsave(plot = integrate_UMAP_condition, 
       filename = file.path(plots_dir,
                            "26_10_2024_integrated_UMAP_condition.png"),
       w=10, h=10, dpi=300)

