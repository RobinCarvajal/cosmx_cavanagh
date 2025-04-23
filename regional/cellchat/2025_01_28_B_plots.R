
# load objects ----------------------------------------------------------------------

# load aldara cellchat object
aldara_cellchat_objects_processed <- readRDS(file.path(data_dir,"regional",
                                                       "aldara_cellchat_objects_processed.RDS"))

# load aldara cellchat object
control_cellchat_objects_processed <- readRDS(file.path(data_dir,"regional",
                                                       "control_cellchat_objects_processed.RDS"))

# creating a path for plots 
plots_dir <- file.path(data_dir,"regional","cellchat", "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Plot 1 ------------------------------------------------------------------------

# aldara
for (region in names(aldara_cellchat_objects_processed)){
  data <- aldara_cellchat_objects_processed[[region]]
  groupSize <- as.numeric(table(data@idents))
  
  # save plot
  region_plots_aldara_dir <- file.path(plots_dir, region, "communication_network")
  dir.create(region_plots_aldara_dir, showWarnings = FALSE, recursive = TRUE)
  
  png(file.path(region_plots_aldara_dir, paste0(region, "_aldara_interaction_plot.png")), width = 4000, height = 2000, res = 300)
  par(mfrow = c(1, 2), xpd = TRUE)  # Set layout for two plots side by side
  netVisual_circle(data@net$count, vertex.weight = rowSums(data@net$count), 
                   weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
  netVisual_circle(data@net$weight, vertex.weight = rowSums(data@net$weight), 
                   weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
  dev.off()  # Save the combined plot
  
}

# control
for (region in names(control_cellchat_objects_processed)){
  data <- control_cellchat_objects_processed[[region]]
  groupSize <- as.numeric(table(data@idents))
  
  # save plot
  region_plots_control_dir <- file.path(plots_dir, region, "communication_network")
  dir.create(region_plots_control_dir, showWarnings = FALSE, recursive = TRUE)
  
  png(file.path(region_plots_control_dir, paste0(region, "_control_interaction_plot.png")), width = 4000, height = 2000, res = 300)
  par(mfrow = c(1, 2), xpd = TRUE)  # Set layout for two plots side by side
  netVisual_circle(data@net$count, vertex.weight = rowSums(data@net$count), 
                   weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
  netVisual_circle(data@net$weight, vertex.weight = rowSums(data@net$weight), 
                   weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
  dev.off()  # Save the combined plot
  
}

# heatmap ------------------------------------------------------------------------

# aldara
for (region in names(aldara_cellchat_objects_processed)){
  data <- aldara_cellchat_objects_processed[[region]]
  
  # save plot
  region_plots_aldara_dir <- file.path(plots_dir, region, "heatmap")
  dir.create(region_plots_aldara_dir, showWarnings = FALSE, recursive = TRUE)
  
  png(file.path(region_plots_aldara_dir, paste0(region, "_aldara_heatmap.png")), width = 2000, height = 2000, res = 300)
  print(netVisual_heatmap(data, measure = "count", color.heatmap = "Blues"))
  dev.off()  # Save the combined plot
  
}

# control 
for (region in names(control_cellchat_objects_processed)){
  data <- control_cellchat_objects_processed[[region]]
  
  # save plot
  region_plots_control_dir <- file.path(plots_dir, region, "heatmap")
  dir.create(region_plots_control_dir, showWarnings = FALSE, recursive = TRUE)
  
  png(file.path(region_plots_control_dir, paste0(region, "_control_heatmap.png")), width = 2000, height = 2000, res = 300)
  print(netVisual_heatmap(data, measure = "count", color.heatmap = "Blues"))
  dev.off()  # Save the combined plot
  
}


