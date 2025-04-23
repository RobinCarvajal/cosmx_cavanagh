

# load brain obj --------------------------------------------------

brain <- readRDS(file.path(data_dir, "objects","brain_objects",
                           "2024_10_26_brain_v1.3_integrated.RDS"))

# load all polygons -------------------------------------------------------

polygons_dir <- file.path(main_dir, "New_Analysis",
                          "single_sample_objects")

single_sample_objects <- list.files(polygons_dir, full.names = TRUE)

# Create a list of loaded objects named by file names
loaded_objects <- lapply(single_sample_objects, function(file) {
  # Replace with appropriate function to load the file (e.g., readRDS, read.csv)
  readRDS(file)  # Example for .RDS files
})

# Name each element in the list with the base file name (without the path and extension)
names(loaded_objects) <- basename(single_sample_objects) %>% tools::file_path_sans_ext()

loaded_objects[["s1055ant_bs1_polygons"]] # example use 

# erase all images from brain obj -----------------------------------------

brain@images <- list()


# change key of all polygons ----------------------------------------------

for (i in 1:length(loaded_objects)){
  
  # Extract the image name and remove "polygons" from the end if it exists
  image_name <- sub("_polygons$", "", names(loaded_objects)[i])
  
  # Extract the polygons
  polygons <- loaded_objects[[i]]
  
  # Change the key
  new_image_key <- paste0(image_name, "_")
  Key(polygons) <- new_image_key
  
  print(Key(polygons))
  
  # Assign the polygons to the brain object
  brain@images[[image_name]] <- polygons
}


# save brain obj ----------------------------------------------------------

saveRDS(brain, file.path(objects_dir, "working_objets",
                         "26_10_2024_brain_v1.4_polygons_corrected.RDS"))
