
# load the brain object 
brain <- readRDS(file.path(objects_dir,"brain_objects",
                           "2024_11_05_brain_v1.5_broad.labels_added.RDS"))


# subsetting brain object -------------------------------------------------

# subset only cell types of interest
celltypes <- c("Astrocytes", "Immune")

# keep cells from the regions of interest
celltypes.sub <- subset(brain, broad.labels %in% celltypes)

# -------------------------------------------------------------------------


# need 1 comparison

# * differences between brain conditions
# Control vs Aldara


# Pseudo-Bulking  ------------------------------------------------

# defining grouping variables to aggregate data
group_vars <- c("broad.labels",
                "mouse_group",
                "slice_position",
                "brain_slice",
                "condition")

# usign custom function
pseudo_brain <- PseudoBulk_fn(obj = celltypes.sub,
                              vars = group_vars)

# save the pseudo object
saveRDS(object = pseudo_brain, 
        file = file.path(data_dir,"broad_single_cell",
                        "broad.labels_pseudo_bulk_condition_astro+immune.RDS"))

