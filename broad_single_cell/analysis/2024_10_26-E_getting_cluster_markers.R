# run in WSL
source("/mnt/c/Users/admin/Documents/GitHub/CosMx_Analysis/Set_up_WSL.R")

# load the integrated object
brain <- readRDS(file.path(objects_dir,"working_objects",
                                     "26_10_2024_brain_v1.3_integrated.RDS"))

# create a dir to store the markers
markers_dir <- file.path(data_dir,"single_cell", "cluster_markers")

# Create the directory if they do not already exist
if (!dir.exists(markers_dir)) {dir.create(markers_dir, recursive = TRUE)}

# get all markers ---------------------------------------------------------

brain.markers <- FindAllMarkers(brain, only.pos = TRUE)
# save the object 
saveRDS(brain.markers,
        file.path(markers_dir,"26_10_2024_all_markers.RDS"))

# save each cluster individually as csv -----------------------------------

for (cluster in unique(brain$integrated_clusters)){
  
  brain.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    ungroup() -> cluster.markers
  
  write.csv(cluster.markers, 
            file.path(markers_dir,
                      paste0("26_10_2024_",cluster,"_markers.csv")))
  
}

