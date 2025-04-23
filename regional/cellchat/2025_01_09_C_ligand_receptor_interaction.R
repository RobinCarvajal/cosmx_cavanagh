
# Load control region.object_list ------------------------------------------
control_region.objects_list <- readRDS(file.path(data_dir,"regional","control_region.objects_list.RDS"))
# Load control region.object_list ------------------------------------------
aldara_region.objects_list <- readRDS(file.path(data_dir,"regional","aldara_region.objects_list.RDS"))

# Defining the database --------------------------------------------
cellchat_DB <- CellChatDB.mouse
# keep only protein interactions
cellchat_DB <- subsetDB(cellchat_DB, non_protein = FALSE) # TRUE when inferring neuron-neuron communication

# Working in parallel ---------------------------------------------------------- 
future::plan("multisession", workers=8)

# Function to create cellchat objects from multiple samples --------------------

createCellChat_multisample <- function(multisample.obj){
  
  #parameter# multisample.obj: list of single sample objects
  #output#                obj: cellchat object
  
  # create data input
  data_input <- createData_multisample(multisample.obj)
  
  # create metadata
  metadata <- createMetadata_multisample(multisample.obj)
  # create levels for fine labels 
  metadata$fine.labels <- as.factor(metadata$fine.labels)
  # create levels for broad labels
  metadata$broad.labels <- as.factor(metadata$broad.labels)
  
  # create spatial locs
  spatial_locs <- createSpatialLocs_multisample(multisample.obj)
  
  # create spatial factors
  spatial_factors <- createSpatialFactors_multisample(metadata, conversion.factor = 1000)
  
  # create cellchat object
  cellchat_obj <- createCellChat(object = data_input, 
                                 datatype = "spatial",
                                 assay="RNA",
                                 meta = metadata, 
                                 group.by = "broad.labels",
                                 coordinates = spatial_locs,
                                 spatial.factors = spatial_factors)
  
  
  # dropping uncessessary levels 
  cellchat_obj <- cellchat_droplevels(cellchat_obj, 
                                      meta_cols = c("broad.labels","fine.labels"), 
                                      clear_idents = TRUE)
  
}

# creating CONTROL/ALDARA objects ----------------------------------------------------

# CONTROL
control_cellchat_objects <- lapply(control_region.objects_list, createCellChat_multisample)
# save the object 
saveRDS(control_cellchat_objects, file = file.path(data_dir, "regional","control_cellchat_objects.RDS"))

# ALDARA
aldara_cellchat_objects <- lapply(aldara_region.objects_list, createCellChat_multisample)
# save the object
saveRDS(aldara_cellchat_objects, file = file.path(data_dir, "regional","aldara_cellchat_objects.RDS"))

# Getting contact range --------------------------------------------------------

# load the brain obj
brain <- readRDS(file.path(objects_dir,"brain_objects","2025_01_22_brain_v1.7_region.labels_added.RDS"))

# get the biggest area for a single cell
brain.meta <- brain@meta.data

# get the max area 
area.median <- median(brain.meta$Area) # in sq microns
# get the diameter
diameter.median <- sqrt(max_area/pi) # in microns
# renaming 
contact_range <- diameter.median

# Getting interaction range ----------------------------------------------------

# I'm going to use 300 microns as the interaction range
# read somewhere that it is how far astrocytes/microglia chemokynes travel 

interaction_range <- 300

# Processing multiple cellchat objects -----------------------------------------

# CONTROL
control_cellchat_objects_processed <- lapply(X = control_cellchat_objects, 
                                             FUN = cellchat_pipeline, 
                                             db=cellchat_DB, 
                                             contact_range= contact_range, 
                                             interaction_range= interaction_range)
# save 
saveRDS(control_cellchat_objects_processed,
        file.path(data_dir,"regional","control_cellchat_objects_processed.RDS"))

# ALDARA
aldara_cellchat_objects_processed <- lapply(X = aldara_cellchat_objects, 
                                            FUN = cellchat_pipeline, 
                                            db=cellchat_DB, 
                                            contact_range= contact_range, 
                                            interaction_range= interaction_range)
# save
saveRDS(aldara_cellchat_objects_processed,
        file.path(data_dir,"regional","aldara_cellchat_objects_processed.RDS"))





