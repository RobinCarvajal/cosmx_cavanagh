
# read barcodes
barcodes <- readRDS(file.path(data_dir,"fine_single_cell",
                            "22_11_2024_fine.labels_barcodes.RDS"))

# read the brain object
brain <- readRDS(file.path(objects_dir,"brain_objects",
                           "05_11_2024_brain_v1.5_broad.labels_added.RDS"))

# get the metadata
brain.meta <- brain@meta.data

# add fine.labels column 
brain.meta$fine.labels <- NA

for (subcluster in names(barcodes)){
  
  sub.barcodes <- barcodes[[subcluster]]
  
  brain.meta$fine.labels[rownames(brain.meta) %in% sub.barcodes] <- subcluster
  
  print(unique(brain.meta$fine.labels))
  
}


# save the new metadata to brain object
brain <- AddMetaData(brain, metadata = brain.meta)


fine.labels_umap <- DimPlot(brain, reduction="umap.integrated",group.by = "fine.labels",
                            label= T, label.box = T)

# save image
ggsave(file.path(data_dir,"fine_single_cell",
                 "22_11_2024_fine.labels_umap.png"),
       fine.labels_umap, width = 14, height = 10, dpi = 600)


# save object
saveRDS(brain, file.path(objects_dir,"brain_objects",
                         "22_11_2024_brain_v1.6_fine.labels_added.RDS"))
