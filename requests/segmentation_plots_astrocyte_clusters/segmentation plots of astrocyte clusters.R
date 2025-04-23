# load brain object
brain <- readRDS(file.path(objects_dir, "brain_objects", 
                            "2025_01_22_brain_v1.7_region.labels_added.RDS"))

# path to save the plots
results_dir <- file.path(disk, "cosmx_results","requests","Lilya")
# crate path 
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

####--- code ---####

# get default colors for the dimplot
library(scales)
cluster_n <- length(unique(brain$integrated_clusters))
dimplot_colors <- hue_pal()(cluster_n)

# cluster colours
cols = dimplot_colors[c(7,11,12)]

# samples
ant_control_sample <- "s1055ant_bs1"
ant_aldara_sample <- "s1055ant_bs3"
pos_control_sample <- "s1055pos_bs1"
pos_aldara_sample <- "s1055pos_bs3"

# amke a subset of brain with only clusters 7,11 and 12
brain_sub <- subset(brain, integrated_clusters %in% c("7","11","12"))
unique(brain_sub$integrated_clusters)

# segmentation plot only with clusters 7, 11 and 12
Idents(brain_sub) <- "integrated_clusters"

# anterior control plot ----------------------------------------------------
plot_ant_control <- ImageDimPlot(brain_sub, fov = ant_control_sample, cols = cols, split.by = "integrated_clusters",
                     boundaries = "segmentation",
                     flip_xy = FALSE,
                     border.size = 0.05)+
  NoLegend()
# save plot as svg
svg(file.path(results_dir, "ant_control_clusters_7_11_12.svg"), width = 20, height = 10)
print(plot_ant_control)
dev.off()


# anterior aldara plot -----------------------------------------------------
plot_ant_aldara <- ImageDimPlot(brain_sub, fov = ant_aldara_sample, cols = cols, split.by = "integrated_clusters",
                                 boundaries = "segmentation",
                                 flip_xy = FALSE,
                                 border.size = 0.05)+
  NoLegend()
# save plot as svg
svg(file.path(results_dir, "ant_aldara_clusters_7_11_12.svg"), width = 20, height = 10)
print(plot_ant_aldara)
dev.off()

# posterior control plot ---------------------------------------------------
plot_pos_control <- ImageDimPlot(brain_sub, fov = pos_control_sample, cols = cols, split.by = "integrated_clusters",
                                  boundaries = "segmentation",
                                  flip_xy = FALSE,
                                  border.size = 0.05)+
  NoLegend()
# save plot as svg
svg(file.path(results_dir, "pos_control_clusters_7_11_12.svg"), width = 20, height = 10)
print(plot_pos_control)
dev.off()

# posterior aldara plot ----------------------------------------------------
plot_pos_aldara <- ImageDimPlot(brain_sub, fov = pos_aldara_sample, cols = cols, split.by = "integrated_clusters",
                                boundaries = "segmentation",
                                flip_xy = FALSE,
                                border.size = 0.05)+
  NoLegend()
# save plot as svg
svg(file.path(results_dir, "pos_aldara_clusters_7_11_12.svg"), width = 20, height = 10)
print(plot_pos_aldara)
dev.off()