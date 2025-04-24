
# load the object 
brain <- readRDS(file.path(data_dir, "objects", "brain_objects",
                           "25_10_2024_brain_v1.0.RDS"))

# quality control 

# set a QC column 
brain$QC <- NA

# set the QC column to Pass for the cells that meet the requirements 

# nFeature_RNA > 0, nCount_RNA > 20, nFeature_RNA < 2000
qc_parameters <- brain$nFeature_RNA > 0 & brain$nCount_RNA > 20 & brain$nFeature_RNA < 2000
brain$QC[qc_parameters] <- "Pass"

# if did not pass set it to Fail
brain$QC[is.na(brain$QC)] <- "Fail"

# how many cells are Pass
table(brain$QC)
# Fail   Pass 
# 58147 755248

# removing cells ----------------------------------------------------------

# keep only cells that are Pass in QC column
brain_qc <- subset(brain, subset = QC == "Pass")

# save the object
saveRDS(brain_qc, file.path(data_dir, "objects", "brain_objects",
                         "25_10_2024_brain_v1.1_QC.RDS"))



# getting cell counts pre and post QC -----------------------------------------

# pre QC

fov_list <- unique(brain$sample_name)

#using dplyr create a table with pass and fail values per sample_name 
counts <- brain@meta.data %>%
  group_by(sample_name) %>%
  summarise(n_cells_preQC = n(),
            n_cells_postQC = n() - sum(QC == "Fail"),
            difference = sum(QC == "Fail"))

# save as csv 
write.csv(counts, file.path(data_dir, "single_cell", "QC_counts.csv"))
