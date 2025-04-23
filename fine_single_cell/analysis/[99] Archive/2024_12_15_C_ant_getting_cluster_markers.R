
# load broad object -------------------------------------------------------

broad <- readRDS(file.path(objects_dir,"broad.label_objects",
                           "2024_12_15_broad_ant_integrated.RDS"))

# create a dir to store the markers
markers_dir <- file.path(data_dir,"fine_single_cell", "ant", "cluster_markers")

# Create the directory if they do not already exist
if (!dir.exists(markers_dir)) {dir.create(markers_dir, recursive = TRUE)}

# names for output files --------------------------------------------------

markers_obj_name <- "2024_12_15_ant_all_broad_markers.RDS"


# process -----------------------------------------------------------------

# empty list where to save the markers
broad.markers <- list()

# getting markers in a loop 

for (cluster in names(broad)){
  
  print(cluster)
  
  cluster.obj <-  broad[[cluster]]
  
  # set idents to the integrated subclusters
  Idents(cluster.obj) <- "integrated_subclusters"
  
  # get the markers
  cluster.obj.markers <- FindAllMarkers(cluster.obj, only.pos = TRUE)
  
  # add the cluster.obj.markers to the broad.markers list 
  broad.markers[[cluster]] <- cluster.obj.markers
  
}

# save the broad.markers object
saveRDS(broad.markers, file.path(markers_dir, markers_obj_name))








