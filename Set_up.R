
# Data Disk ----------------------------------------------------------

disk <- "F:"

# Scripts Directory -------------------------------------------------------

scripts_dir <- file.path(disk, "cosmx_scripts")


# Set up Directory --------------------------------------------------------

setup_dir <- file.path(scripts_dir,"set_up")


# Data Directories --------------------------------------------------------

data_dir <- file.path(disk, "cosmx_data")

##### objects #####
objects_dir <- file.path(data_dir, "objects")

##### single cell #####
single_cell_dir <- file.path(data_dir, "single_cell")

##### niches #####
niches_dir <- file.path(data_dir, "niches")



#### running scripts 

# Functions
source(file.path(setup_dir, "Functions.R"))
source(file.path(setup_dir, "FN_differential_expression.R"))
source(file.path(setup_dir, "FN_cellchat.R"))
source(file.path(setup_dir, "SL2_Functions.R"))

# Libraries
source(file.path(setup_dir, "Libraries.R"))




