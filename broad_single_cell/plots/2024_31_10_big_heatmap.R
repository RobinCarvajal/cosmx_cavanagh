
# load the integrated object
brain <- readRDS(file.path(data_dir,"objects","working_objects",
                           "26_10_2024_brain_v1.3_integrated.RDS"))

# load markers
brain.markers <- readRDS(file.path(data_dir,"single_cell","cluster_markers",
                                   "26_10_2024_all_markers.RDS"))

# create plots dir --------------------------------------------------------

# create a dir to store the plots
plots_dir <- file.path(data_dir,"single_cell", "plots")

# Create the directory if they do not already exist
if (!dir.exists(plots_dir)) {dir.create(plots_dir, recursive = TRUE)}

# Heatmap of top 10 marker genes 
brain.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.58) %>% # 1.5 times fold change
  slice_head(n = 10) %>%
  ungroup() -> top10

# need to downsample brain obj (faster)
brain_downsampled <- subset(x = brain, downsample = 100)

# plot
top10_heatmap <- DoHeatmap(brain_downsampled, features = top10$gene) + NoLegend()

ggsave(plot = top10_heatmap,
       filename = file.path(plots_dir,"31_10_20204_top_10_marker_genes_heatmap.png"),
       width = 10, height = 15, dpi = 600)
