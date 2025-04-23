
brain.list <- readRDS(file.path(data_dir,"niches","2024_11_28_brain.list_v1.0.RDS"))


# trying function in a loop  ----------------------------------------------


# all samples
all_samples <- names(brain.list)
# anterior samples
anterior_samples <- all_samples[grep("ant", all_samples)]
# posterior samples
posterior_samples <- all_samples[grep("pos", all_samples)]
  

# niches list -------------------------------------------------------------

# create empty list 
ant_niches.list <- list()
pos_niches.list <- list()

# ANTERIOR samples (5 niches) --------------------------------------------------------

# processing 

for (sample in anterior_samples) {
  sample_obj <- brain.list[[sample]]
  
  # Define the range of niches.k values
  niches_k_values <- 5:9
  
  # Create assays for each niches.k value
  for (niches_k in niches_k_values) {
    
    cluster_name <- paste0("niches_k", niches_k)
    
    sample_obj <- BuildNicheAssay(
      object = sample_obj, 
      fov = sample, 
      group.by = "fine.labels",
      niches.k = niches_k, 
      neighbors.k = 30, 
      cluster.name = cluster_name
    )
  }
  
  # Store the updated sample object in the list
  ant_niches.list[[sample]] <- sample_obj
}

# save
saveRDS(ant_niches.list, file.path(data_dir,"niches","2024_12_04_niches.list_ant_k_5_9.RDS"))


# POSTERIOR samples (5:7 niches) --------------------------------------------------------

pos_niches.list <- list()
# processing 
for (sample in posterior_samples) {
  sample_obj <- brain.list[[sample]]
  
  # Define the range of niches.k values
  niches_k_values <- 5:7
  
  # Create assays for each niches.k value
  for (niches_k in niches_k_values) {
    
    cluster_name <- paste0("niches_k", niches_k)
    
    sample_obj <- BuildNicheAssay(
      object = sample_obj, 
      fov = sample, 
      group.by = "fine.labels",
      niches.k = niches_k, 
      neighbors.k = 30, 
      cluster.name = cluster_name
    )
  }
  
  # Store the updated sample object in the list
  pos_niches.list[[sample]] <- sample_obj
}

# save 

saveRDS(pos_niches.list, file.path(data_dir,"niches","2024_12_17_niches.list_pos_k_5_7.RDS"))


# POSTERIOR samples (8,12 niches) --------------------------------------------------------

# processing 

for (sample in posterior_samples) {
  sample_obj <- brain.list[[sample]]
  
  # Define the range of niches.k values
  niches_k_values <- 8:12
  
  # Create assays for each niches.k value
  for (niches_k in niches_k_values) {
    
    cluster_name <- paste0("niches_k", niches_k)
    
    sample_obj <- BuildNicheAssay(
      object = sample_obj, 
      fov = sample, 
      group.by = "fine.labels",
      niches.k = niches_k, 
      neighbors.k = 30, 
      cluster.name = cluster_name
    )
  }
  
  # Store the updated sample object in the list
  pos_niches.list[[sample]] <- sample_obj
}

# save 

saveRDS(pos_niches.list, file.path(data_dir,"niches","2024_12_04_niches.list_pos_k_8_12.RDS"))
