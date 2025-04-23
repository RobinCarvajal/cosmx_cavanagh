# load the brain obj with regions annotated
brain <- readRDS(file.path(objects_dir,"brain_objects",
                         "2025_01_22_brain_v1.7_region.labels_added.RDS"))


# subset only needed regions 
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

# add cells from anterior caudoputamen - reason: overlap with 
ant_caudoputamen.cells <- brain@meta.data%>%
  filter(regions.major == "caudoputamen" & slice_position == "ant") %>%
  row.names()

# update roi cells 
roi.cells <- c(roi.cells, ant_caudoputamen.cells)
  
# subset the object with only roi cells 
roi.sub <- subset(brain, cells = roi.cells)




# pseudo bulking 
group_vars <- c("regions.major",
                "mouse_group",
                "slice_position",
                "brain_slice",
                "condition")

regions_bulk <- PseudoBulk_fn(roi.sub, vars = group_vars)

# save the pseudo bulk object
saveRDS(regions_bulk, file.path(data_dir,"regional",
                                "2025_01_23_region-bulk.RDS"))