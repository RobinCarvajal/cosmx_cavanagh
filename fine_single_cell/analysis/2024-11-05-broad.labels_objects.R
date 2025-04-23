
# load object -------------------------------------------------------------

brain <- readRDS(file.path(objects_dir,"brain_objects",
                           "05_11_2024_brain_v1.5_broad.labels_added.RDS"))


# make subsets per broad label  ----------------------------------------------

# get the broad labels
broad_labels <- unique(brain$broad.labels)

# create a list of subsets
broad <- list()
for (label in broad_labels){
  subset <- subset(brain, broad.labels == label)
  #add to the list
  broad[[label]] <- subset
  print(label)
}

# save object ----------------------------------------------------------

saveRDS(broad, file.path(objects_dir,"broad.label_objects",
                          "05_11_2024_broad_v1.0.RDS"))




















