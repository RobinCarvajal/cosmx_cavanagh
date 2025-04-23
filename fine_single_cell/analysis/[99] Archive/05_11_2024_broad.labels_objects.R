
# load object -------------------------------------------------------------

brain <- readRDS(file.path(objects_dir,"brain_objects",
                           "2024_11_05_brain_v1.5_broad.labels_added.RDS"))


# separate by position  ----------------------------------------------------
all_samples <- unique(brain$sample_name)
anterior_samples <- grep("ant", all_samples, value=T)
posterior_samples <- grep("pos", all_samples, value=T)

# subset the data based on position of samples
brain_ant <- subset(brain, sample_name %in% anterior_samples)
brain_pos <- subset(brain, sample_name %in% posterior_samples)

# get all the broad.labels ------------------------------------------------

broad.labels <- unique(brain$broad.labels)
# celltypes excluded du to small size
celltypes_to_exclude <- c("cluster0","cluster1","Endothelial","Ependymal")

# filtering celltypes
filtered.labels <- broad.labels[!broad.labels %in% celltypes_to_exclude]
# grop levels not present
filtered.labels <- droplevels(filtered.labels)

# (ANTERIOR) make subsets per broad label  -----------------------------------

# create a list of subsets
broad_ant <- list()
for (label in filtered.labels){
  subset <- subset(brain_ant, broad.labels == label)
  #add to the list
  broad_ant[[label]] <- subset
  print(label)
}

# save object
saveRDS(broad_ant, file.path(objects_dir,"broad.label_objects",
                         "2024_12_15_broad_ant.RDS"))


# (POSTERIOR) make subsets per broad label  -----------------------------------

# create a list of subsets
broad_pos <- list()
for (label in filtered.labels){
  subset <- subset(brain_pos, broad.labels == label)
  #add to the list
  broad_pos[[label]] <- subset
  print(label)
}

# save object
saveRDS(broad_pos, file.path(objects_dir,"broad.label_objects",
                         "2024_12_15_broad_pos.RDS"))

















