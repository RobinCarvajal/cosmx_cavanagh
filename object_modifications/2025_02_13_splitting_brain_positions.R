
# load data
brain <- readRDS(file.path(data_dir,"objects", "brain_objects", 
                           "2025_01_22_brain_v1.7_region.labels_added.RDS"))

# Create two objects - Split by position 
brain_anterior <- subset(brain, subset = slice_position == "ant")
brain_posterior <- subset(brain, subset = slice_position == "pos")

# save the RDS objects
saveRDS(brain_anterior, file.path(data_dir,"objects", "brain_objects", 
                                  "2025_02_13_brain_v1.7_anterior.RDS"))
saveRDS(brain_posterior, file.path(data_dir,"objects", "brain_objects", 
                                   "2025_02_13_brain_v1.7_posterior.RDS"))