
# Load the brain object with annotated regions
brain <- readRDS(file.path(objects_dir,"brain_objects","2025_01_22_brain_v1.7_region.labels_added.RDS"))

# subset the brain to include only Regions of Interest ------------------------

# regions of interest
roi <- c("anterior_cingulate_area",
         #"caudoputamen",
         "nucleus_accumbens",
         "thalamus",
         "hippocampus",
         "motor_areas",
         "agranular_insular_area",
         "amygdala")

# keep cells from the regions of interest
roi.sub <- subset(brain, regions.major %in% roi)

roi.cells <- brain@meta.data %>%
  filter(regions.major %in% roi) %>%
  row.names()

# add cells from anterior caudoputamen
ant_caudoputamen.cells <- brain@meta.data%>%
  filter(regions.major == "caudoputamen" & slice_position == "ant") %>%
  row.names()

# update roi cells 
roi.cells <- c(roi.cells, ant_caudoputamen.cells)

# subset the object with only roi cells 
roi.sub <- subset(brain, cells = roi.cells)


# other stuff ------------------
# list of aldara samples 
control_samples <- unique(subset(brain@meta.data, condition=="C")$sample_name)
# list of control samples 
aldara_samples <- unique(subset(brain@meta.data, condition=="A")$sample_name)

# Control brain subset
control_brain <- subset(roi.sub, sample_name %in% control_samples)
# Aldata brain subset 
aldara_brain <- subset(roi.sub, sample_name %in% aldara_samples)

# custom function --------------------------------------------------------------

create_region.object_list <- function(obj){
  
  # empty list to store all the regions with multiple samples
  region.objects_list <- list()
  
  # processing loop
  for (region in unique(obj$regions.major)){
    
    # skip "not assigned"
    if (region == "not_assigned") next
    
    # sample subset
    region.sub <- subset(obj, regions.major == region)
    
    # empty list to store multiple sample objects per region
    samples <- list()
    for (sample in unique(region.sub$sample_name)){
      
      sample.sub <- subset(region.sub, sample_name == sample)
      # adding the region subset to the corresponding sample
      samples[[sample]] <- sample.sub
    }
    
    # save the sample objects to the corresponding region
    region.objects_list[[region]] <- samples
    
  }
  
  return(region.objects_list)
  
}

# CONTROL samples --------------------------------------------------------------

# using custom function 
control_region.objects_list <- create_region.object_list(control_brain)
# saving object
saveRDS(control_region.objects_list, file.path(data_dir,"regional","control_region.objects_list.RDS"))

# loop for ALDARA samples --------------------------------------------------------------

# using custom function 
aldara_region.objects_list <- create_region.object_list(aldara_brain)
# saving object
saveRDS(aldara_region.objects_list, file.path(data_dir,"regional","aldara_region.objects_list.RDS"))
