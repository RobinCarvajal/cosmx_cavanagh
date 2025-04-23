

# load the slide objects --------------------------------------------------

s1055ant <- readRDS(file.path(main_dir, "Objects_NO_QC", "original_1055_1.RDS"))
s1055pos <- readRDS(file.path(main_dir, "Objects_NO_QC", "original_1055_2.RDS"))
s1056ant <- readRDS(file.path(main_dir, "Objects_NO_QC", "original_1056_1.RDS"))
s1056pos <- readRDS(file.path(main_dir, "Objects_NO_QC", "original_1056_2.RDS"))

# needed to process data 
library(future)
options(future.globals.maxSize = 8 * 1024^3)  # Set to 8 GiB


# change cell names -------------------------------------------------------

colnames(s1055ant) <- s1055ant$brain_cell_id
colnames(s1055pos) <- s1055pos$brain_cell_id
colnames(s1056ant) <- s1056ant$brain_cell_id
colnames(s1056pos) <- s1056pos$brain_cell_id


# merging objects ---------------------------------------------------------

brain <- merge(s1055ant, y=c(s1055pos, s1056ant, s1056pos))


# renaming the images key -------------------------------------------------

Key(brain@images[["Glas1055_1"]]) <- "s1055ant_"
Key(brain@images[["Glas1055_2"]]) <- "s1055pos_"
Key(brain@images[["Glas1056_1"]]) <- "s1056ant_"
Key(brain@images[["Glas1056_2"]]) <- "s1056pos_"

# rename the images
names(brain@images) <- c("s1055ant", "s1055pos", "s1056ant", "s1056pos")


# adding condition --------------------------------------------------------


# if bs1 or bs2, then condition is control(C), otherwise is Aldara(A)
brain$condition <- ifelse(grepl("bs1|bs2", brain$brain_slice), "C", "A")

# setting levels
brain$condition <- factor(brain$condition, levels = c("C","A"))

# check levels (control first, aldara second)
levels(brain$condition)

# checking
table(brain$condition,brain$brain_slice)
#      bs1    bs2    bs3    bs4
# A      0      0 191360 188526
# C 215103 218406      0      0


# save the object ---------------------------------------------------------
saveRDS(brain, file = file.path(main_dir, "New_Analysis","objects",
                                "25_10_2024_brain_v1.0.RDS"))

