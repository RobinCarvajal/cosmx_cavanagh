
# Glas1055_1 --------------------------------------------------------------

# load object
org_1055_1 <- readRDS("F:/CosMx_Analysis/Flat_Files/Slide_1055_1/Glas1055_1/slide_objects/Glas1055_1_polygons_seuratObject.RDS")

##### defining mouse_group #####
org_1055_1$mouse_group <- "s1055"

##### defining slice_position #####
org_1055_1$slice_position <- "ant"

##### defining brain_slice #####

# Define the indices for each brain slice
brainslice_1 <- 1:73      # TL
brainslice_2 <- 74:169    # TR
brainslice_3 <- c(170:249, 318, 319)  # BL
brainslice_4 <- 250:317    # BR

# Create a vector to hold the brain slice labels
brain_slices_1055_1 <- rep(NA, length(org_1055_1$cell_id))

# Assign labels based on the fov values in the metadata
brain_slices_1055_1[org_1055_1@meta.data$fov %in% brainslice_1] <- "bs1"
brain_slices_1055_1[org_1055_1@meta.data$fov %in% brainslice_2] <- "bs2"
brain_slices_1055_1[org_1055_1@meta.data$fov %in% brainslice_3] <- "bs3"
brain_slices_1055_1[org_1055_1@meta.data$fov %in% brainslice_4] <- "bs4"

# Add the brain slice information to the metadata of the Seurat object
org_1055_1$brain_slice <- brain_slices_1055_1

# check original cell counts
sum(org_1055_1$brain_slice == "bs1") # 37599
sum(org_1055_1$brain_slice == "bs2") # 46855
sum(org_1055_1$brain_slice == "bs3") # 42269
sum(org_1055_1$brain_slice == "bs4") # 35761

# Check for NA values in the brain_slice column
na_count <- sum(is.na(org_1055_1$brain_slice))

##### defining sample_name #####
# defining by merging mouse_group, slice_position and brain_slice

# First merge mouse_group and slice_position without a separator
temp_name <- paste(org_1055_1$mouse_group, org_1055_1$slice_position, sep = "")

# Then merge the result with brain_slice using "_"
org_1055_1$sample_name <- paste(temp_name, org_1055_1$brain_slice, sep = "_")

##### defining brain_cell_id #####
# defining by merging sample_name and cell_id
org_1055_1$brain_cell_id <- paste(org_1055_1$sample_name, org_1055_1$cell_id, sep = "_")

##### save RDS #####
saveRDS(org_1055_1, file.path(main_dir,"Objects_NO_QC","original_1055_1.RDS"))

# Glas1055_2 --------------------------------------------------------------

# load object
org_1055_2 <- readRDS("F:/CosMx_Analysis/Flat_Files/Slide_1055_2/Glas1055_2/Glas1055_2_polygons_seuratObject.RDS")

##### defining mouse_group #####
org_1055_2$mouse_group <- "s1055"

##### defining slice_position #####
org_1055_2$slice_position <- "pos"

##### defining brain_slice #####

# Define the indices for each brain slice
brainslice_1 <- 1:92      # TL
brainslice_2 <- c(93:202,382:386)    # TR
brainslice_3 <- c(203:290,387:390)  # BL
brainslice_4 <- c(291:381)    # BR

# Create a vector to hold the brain slice labels
brain_slices_1055_2 <- rep(NA, length(org_1055_2$cell_id))

# Assign labels based on the fov values in the metadata
brain_slices_1055_2[org_1055_2@meta.data$fov %in% brainslice_1] <- "bs1"
brain_slices_1055_2[org_1055_2@meta.data$fov %in% brainslice_2] <- "bs2"
brain_slices_1055_2[org_1055_2@meta.data$fov %in% brainslice_3] <- "bs3"
brain_slices_1055_2[org_1055_2@meta.data$fov %in% brainslice_4] <- "bs4"

# Add the brain slice information to the metadata of the Seurat object
org_1055_2$brain_slice <- brain_slices_1055_2

# check original cell counts
sum(org_1055_2$brain_slice == "bs1", na.rm = T) # 53197
sum(org_1055_2$brain_slice == "bs2", na.rm = T)  # 68856
sum(org_1055_2$brain_slice == "bs3", na.rm = T)  # 59008
sum(org_1055_2$brain_slice == "bs4", na.rm = T)  # 59497

# Check for NA values in the brain_slice column
na_count <- sum(is.na(org_1055_2$brain_slice))
print(na_count)

##### defining sample_name #####
# defining by merging mouse_group, slice_position and brain_slice

# First merge mouse_group and slice_position without a separator
temp_name <- paste(org_1055_2$mouse_group, org_1055_2$slice_position, sep = "")

# Then merge the result with brain_slice using "_"
org_1055_2$sample_name <- paste(temp_name, org_1055_2$brain_slice, sep = "_")

##### defining brain_cell_id #####
# defining by merging sample_name and cell_id
org_1055_2$brain_cell_id <- paste(org_1055_2$sample_name, org_1055_2$cell_id, sep = "_")

##### save RDS #####
saveRDS(org_1055_2, file.path(main_dir,"Objects_NO_QC","original_1055_2.RDS"))

# Glas1056_1 --------------------------------------------------------------

# load object
org_1056_1 <- readRDS("F:/CosMx_Analysis/Flat_Files/Slide_1056_1/Glas1056_1/Glas1056_1_polygons_seuratObject.RDS")

##### defining mouse_group #####
org_1056_1$mouse_group <- "s1056"

##### defining slice_position #####
org_1056_1$slice_position <- "ant"

##### defining brain_slice #####

# Define the indices for each brain slice
brainslice_1 <- 1:82      # TL
brainslice_2 <- 83:159    # TR
brainslice_3 <- 160:227  # BL
brainslice_4 <- 228:305    # BR

# Create a vector to hold the brain slice labels
brain_slices_1056_1 <- rep(NA, length(org_1056_1$cell_id))

# Assign labels based on the fov values in the metadata
brain_slices_1056_1[org_1056_1@meta.data$fov %in% brainslice_1] <- "bs1"
brain_slices_1056_1[org_1056_1@meta.data$fov %in% brainslice_2] <- "bs2"
brain_slices_1056_1[org_1056_1@meta.data$fov %in% brainslice_3] <- "bs3"
brain_slices_1056_1[org_1056_1@meta.data$fov %in% brainslice_4] <- "bs4"

# Add the brain slice information to the metadata of the Seurat object
org_1056_1$brain_slice <- brain_slices_1056_1

# check original cell counts
sum(org_1056_1$brain_slice == "bs1") # 56459
sum(org_1056_1$brain_slice == "bs2") # 45241
sum(org_1056_1$brain_slice == "bs3") # 43433
sum(org_1056_1$brain_slice == "bs4") # 43342

# Check for NA values in the brain_slice column
na_count <- sum(is.na(org_1056_1$brain_slice))
print(na_count)

##### defining sample_name #####
# defining by merging mouse_group, slice_position and brain_slice

# First merge mouse_group and slice_position without a separator
temp_name <- paste(org_1056_1$mouse_group, org_1056_1$slice_position, sep = "")

# Then merge the result with brain_slice using "_"
org_1056_1$sample_name <- paste(temp_name, org_1056_1$brain_slice, sep = "_")

##### defining brain_cell_id #####
# defining by merging sample_name and cell_id
org_1056_1$brain_cell_id <- paste(org_1056_1$sample_name, org_1056_1$cell_id, sep = "_")

##### save RDS #####
saveRDS(org_1056_1, file.path(main_dir,"Objects_NO_QC","original_1056_1.RDS"))


# Glas1056_2 --------------------------------------------------------------

# load object
org_1056_2 <- readRDS("F:/CosMx_Analysis/Flat_Files/Slide_1056_2/Glas1056_2/Glas1056_2_polygons_seuratObject.RDS")

##### defining mouse_group #####
org_1056_2$mouse_group <- "s1056"

##### defining slice_position #####
org_1056_2$slice_position <- "pos"

##### defining brain_slice #####

# Define the indices for each brain slice
brainslice_1 <- c(1:84,338:342)      # TL
brainslice_2 <- c(85:177)    # TR
brainslice_3 <- c(178:246,332:337)  # BL
brainslice_4 <- 247:331    # BR

# Create a vector to hold the brain slice labels
brain_slices_1056_2 <- rep(NA, length(org_1056_2$cell_id))

# Assign labels based on the fov values in the metadata
brain_slices_1056_2[org_1056_2@meta.data$fov %in% brainslice_1] <- "bs1"
brain_slices_1056_2[org_1056_2@meta.data$fov %in% brainslice_2] <- "bs2"
brain_slices_1056_2[org_1056_2@meta.data$fov %in% brainslice_3] <- "bs3"
brain_slices_1056_2[org_1056_2@meta.data$fov %in% brainslice_4] <- "bs4"

# Add the brain slice information to the metadata of the Seurat object
org_1056_2$brain_slice <- brain_slices_1056_2

# check original cell counts
sum(org_1056_2$brain_slice == "bs1") # 67848
sum(org_1056_2$brain_slice == "bs2") # 57454
sum(org_1056_2$brain_slice == "bs3") # 46650
sum(org_1056_2$brain_slice == "bs4") # 49926

# Check for NA values in the brain_slice column
na_count <- sum(is.na(org_1056_2$brain_slice))
print(na_count)

##### defining sample_name #####
# defining by merging mouse_group, slice_position and brain_slice

# First merge mouse_group and slice_position without a separator
temp_name <- paste(org_1056_2$mouse_group, org_1056_2$slice_position, sep = "")

# Then merge the result with brain_slice using "_"
org_1056_2$sample_name <- paste(temp_name, org_1056_2$brain_slice, sep = "_")

##### defining brain_cell_id #####
# defining by merging sample_name and cell_id
org_1056_2$brain_cell_id <- paste(org_1056_2$sample_name, org_1056_2$cell_id, sep = "_")

##### save RDS #####
saveRDS(org_1056_2, file.path(main_dir,"Objects_NO_QC","original_1056_2.RDS"))










             
             
             