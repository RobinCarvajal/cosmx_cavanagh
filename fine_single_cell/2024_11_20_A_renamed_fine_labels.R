
#  load object ------------------------------------------------------------

broad <- readRDS(file.path(objects_dir,"broad.label_objects",
                           "2024_11_06_broad_v1.3_integrated.RDS"))

# mapping labels to clusters ----------------------------------------------

# add a new column to the metadata 

for (celltype in names(broad)){
  

  sub <- broad[[celltype]]
  
  sub.meta<- sub@meta.data
  
  # add the fine labels column
  if (celltype == "Ependymal"){
    sub.meta$fine.labels <- sub.meta$broad.labels
  } else {
    sub.meta$fine.labels <- paste0(sub.meta$broad.labels,"_", sub.meta$integrated_subclusters)
  }
  
  # save the new metadata to sub object
  sub <- AddMetaData(sub, metadata = sub.meta)
  
  print(celltype)
  print(head(sub@meta.data,5))
  
  # add back to broad object
  broad[[celltype]] <- sub
  
  
}


#  get a list of barcodes per fine label ----------------------------------

# barcodes list
barcodes <- list()

for (cluster in names(broad)){
  
  cluster.obj <- broad[[cluster]]
  
  print(cluster)
  
  # list of subclusters
  subclusters <- unique(cluster.obj$fine.labels)
  #print(subclusters)
  
  for (subcluster in subclusters){
    
    # get a subcluster obj
    subcluster.obj <- subset(cluster.obj, fine.labels == subcluster)
    
    
    #print(subcluster.obj)
    # adding to the barcodes list
    barcodes[[subcluster]] <- Cells(subcluster.obj)
    
  }
  
}  

# save the obj ------------------------------------------------------------

# create dir 
dir.create(file.path(data_dir,"fine_single_cell"), showWarnings = FALSE)

saveRDS(barcodes, file.path(data_dir,"fine_single_cell",
                         "22_11_2024_fine.labels_barcodes.RDS"))


