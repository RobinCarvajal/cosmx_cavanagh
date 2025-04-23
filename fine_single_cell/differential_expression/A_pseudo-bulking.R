
# load the brain object 
brain <- readRDS(file.path(objects_dir,"brain_objects",
                           "2024_11_22_brain_v1.6_fine.labels_added.RDS"))

# need 1 comparison

# * differences between brain conditions
# Control vs Aldara


# Pseudo-Bulking  ------------------------------------------------

# defining grouping variables to aggregate data
group_vars <- c("fine.labels",
                "mouse_group",
                "slice_position",
                "brain_slice",
                "condition")

# usign custom function
pseudo_brain <- PseudoBulk_fn(obj = brain,
                              vars = group_vars)

# save the pseudo object
saveRDS(object = pseudo_brain, 
        file = file.path(data_dir,"de_utils",
                        "fine.labels_pseudo_bulk.RDS"))

