
# load broad object -------------------------------------------------------

broad <- readRDS(file.path(objects_dir,"broad.label_objects",
                           "06_11_2024_broad_v1.3_integrated.RDS"))

# create a dir to store the markers
markers_dir <- file.path(data_dir,"broad_single_cell", "cluster_markers")

# Create the directory if they do not already exist
if (!dir.exists(markers_dir)) {dir.create(markers_dir, recursive = TRUE)}


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
saveRDS(broad.markers, file.path(markers_dir,
                                 "06_11_2024_all_broad_markers.RDS"))


system("shutdown /s /f /t 300")






