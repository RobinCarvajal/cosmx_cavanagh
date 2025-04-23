
# Load niches list --------------------------------------------------------
ant_niches.list <- readRDS(file.path(data_dir,"niches","2024_12_04_niches.list_ant_k_5_9.RDS"))
# 8-12
pos_niches.list <- readRDS(file.path(data_dir,"niches","2024_12_04_niches.list_pos_k_8_12.RDS"))
# 5-7
pos_niches.list <- readRDS(file.path(data_dir,"niches","2024_12_04_niches.list_pos_k_5_7.RDS"))

# sample names ------------------------------------------------------------------

# anterior samples
anterior_samples <- names(ant_niches.list)
# posterior samples
posterior_samples <- names(pos_niches.list)

# create dirs  ------------------------------------------------------------

# create folders for each sample 
for(sample in names(niches.list)){
  dir.create(file.path(data_dir,"niches",sample), showWarnings = F)
}


# ANTERIOR niches --------------------------------------------------------

for (sample in anterior_samples){
  
  dir.create(file.path(data_dir,"niches",sample), showWarnings = F)
  
  sample_obj <- ant_niches.list[[sample]]
  
  niches_k <- paste0("niches_k",5)
  
  # set idents as niches
  Idents(sample_obj) <- niches_k
  
  ## all niches plot 
  niche.plot <- ImageDimPlot(sample_obj, 
                             fov=sample,
                             boundaries="segmentation",
                             group.by = niches_k, 
                             border.size = 0.05,
                             border.color = "white",
                             dark.background = T,
                             flip_xy = F) +
    ggtitle(sample)
  
  # save
  ggsave(file.path(data_dir,"niches",sample, paste0(sample,"_niches.png")), 
         plot = niche.plot, width = 10, height = 10, dpi=300)
  
  ## split plot
  niches_split.plot <- ImageDimPlot(sample_obj, 
                                    fov=sample,
                                    boundaries="centroids",
                                    split.by = niches_k, 
                                    #border.size = 0.05,
                                    dark.background = T,
                                    flip_xy = F,
                                    combine = T ) +
    ggtitle(sample)
  
  # save
  ggsave(file.path(data_dir,"niches", sample, paste0(sample,"_all_niches.png")), 
         plot = niches_split.plot, width = 10, height = 10, dpi=300)
  
  
  
}


# POSTERIOR niches --------------------------------------------------------
for (sample in posterior_samples){
  
  dir.create(file.path(data_dir,"niches",sample), showWarnings = F)
  
  sample_obj <- pos_niches.list[[sample]]
  
  niches_k <- paste0("niches_k",7) # <- change the number depending on no. clusters wanted
  
  # set idents as niches
  Idents(sample_obj) <- niches_k
  
  ## all niches plot 
  niche.plot <- ImageDimPlot(sample_obj, 
                             fov=sample,
                             boundaries="segmentation",
                             group.by = niches_k, 
                             border.size = 0.05,
                             border.color = "white",
                             dark.background = T,
                             flip_xy = F) +
    ggtitle(sample)
  
  # save
  ggsave(file.path(data_dir,"niches",sample, paste0(sample,"_niches.png")), 
         plot = niche.plot, width = 10, height = 10, dpi=300)
  
  ## split plot
  niches_split.plot <- ImageDimPlot(sample_obj, 
                                    fov=sample,
                                    boundaries="centroids",
                                    split.by = niches_k, 
                                    #border.size = 0.05,
                                    dark.background = T,
                                    flip_xy = F,
                                    combine = T ) +
    ggtitle(sample)
  
  # save
  ggsave(file.path(data_dir,"niches", sample, paste0(sample,"_all_niches.png")), 
         plot = niches_split.plot, width = 10, height = 10, dpi=300)
  
}