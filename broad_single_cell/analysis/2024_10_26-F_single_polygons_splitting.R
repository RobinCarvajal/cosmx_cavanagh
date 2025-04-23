# dirs
new_analysis_dir <- file.path(main_dir, "New_Analysis")
objects_dir <- file.path(new_analysis_dir, "objects")
single_sample_obj_dir <- file.path(new_analysis_dir, "single_sample_objects")


# load the object 
# using QC object because its lighter to load
brain <- readRDS(file.path(objects_dir,"26_10_2024_brain_v1.3_integrated.RDS")) 

# get all sample names
sample_names <- unique(brain$sample_name)

# for loop to split and save the objects
 for (sample in sample_names){
   
   # make a subset 
   single_sample <- subset(brain, subset = sample_name == sample)
   
   # extract polygons
   image_name <- paste0(unique(single_sample$mouse_group),
                        unique(single_sample$slice_position))
   single_sample_polygons <- single_sample@images[[image_name]]
   
   # changing the Key of the polygons object
   new_image_key <- paste0(image_name,"_",unique(single_sample$brain_slice),"_")
   Key(single_sample@images[[image_name]]) <- new_image_key
   
   print(Key(single_sample@images[[image_name]]))
   
   # save the single sample polygons obj
   saveRDS(single_sample_polygons, file.path(single_sample_obj_dir, paste0(sample,"_polygons.RDS")))
   
   print(sprintf("sample %s saved", sample))
   
 }




















