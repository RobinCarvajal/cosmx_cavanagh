
# load the brain object ---------------------------------------------------

brain <- readRDS(file.path(objects_dir,"brain_objects",
                           "2024_11_22_brain_v1.6_fine.labels_added.RDS"))


# Creating the brain.list object ------------------------------------------

# empty list to save sample objects
brain.list <- list()

# samples vector
samples <- unique(brain$sample_name)

for (sample in samples){
  
  # creatign sample obj
  sample_obj <- subset(brain, sample_name == sample)
  
  # add to the list
  brain.list[[sample]] <- sample_obj
  
}

# save the object 
saveRDS(brain.list, file.path(data_dir,"niches",
                              "2024_11_28_brain.list_v1.0.RDS"))
