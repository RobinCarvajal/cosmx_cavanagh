
# extracting interaction tables

# Aldara ----------------------------------------------------------------------

# load aldara cellchat object
aldara_cellchat_objects_processed <- readRDS(file.path(data_dir,"regional",
                                                       "aldara_cellchat_objects_processed.RDS"))

# Extracting interaction tables
aldara_interaction_tables <- list()
# loop
for (region in names(aldara_cellchat_objects_processed)){
  data <- aldara_cellchat_objects_processed[[region]]
  interactions <- subsetCommunication(data)
  aldara_interaction_tables[[region]] <- interactions
}


# create a interactions folder
interactions_aldara_dir <- file.path(data_dir,"regional","cellchat", "interactions","aldara")
dir.create(interactions_aldara_dir, showWarnings = FALSE, recursive = TRUE)

# save each table as a tsv file 
for (region in names(aldara_interaction_tables)){
  data <- aldara_interaction_tables[[region]]
  write.table(data, file = file.path(interactions_aldara_dir, paste0(region,"_interaction_table.tsv")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# CONTROL ----------------------------------------------------------------------

# load control cellchat object
control_cellchat_objects_processed <- readRDS(file.path(data_dir,"regional",
                                                        "control_cellchat_objects_processed.RDS"))

# Extracting interaction tables
control_interaction_tables <- list()
# loop
for (region in names(control_cellchat_objects_processed)){
  data <- control_cellchat_objects_processed[[region]]
  interactions <- subsetCommunication(data)
  control_interaction_tables[[region]] <- interactions
}

# create a interactions folder
interactions_control_dir <- file.path(data_dir,"regional","cellchat", "interactions","control")
dir.create(interactions_control_dir, showWarnings = FALSE, recursive = TRUE)

# save each table as a tsv file
for (region in names(control_interaction_tables)){
  data <- control_interaction_tables[[region]]
  write.table(data, file = file.path(interactions_control_dir, paste0(region,"_interaction_table.tsv")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
}

