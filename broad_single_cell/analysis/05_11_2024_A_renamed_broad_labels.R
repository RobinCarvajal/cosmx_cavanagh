
#  load object ------------------------------------------------------------

brain <- readRDS(file.path(objects_dir,"brain_objects",
                         "26_10_2024_brain_v1.4_polygons_corrected.RDS"))


# mapping labels to clusters ----------------------------------------------

# set integrated_clusters labels as idents
Idents(brain) <- "integrated_clusters"

# rename the Idents

brain <- RenameIdents(brain, 
                      "0"="cluster0",
                      "1"="cluster1",
                      "2"="Glutamatergic",
                      "3"="GABAergic",
                      "4"="Neurovascular",
                      "5"="Oligodendrocytes",
                      "6"="Glutamatergic",
                      "7"="Astrocytes",
                      "8"="MSNs",
                      "9"="Oligodendrocytes",
                      "10"="Immune",
                      "11"="Astrocytes",
                      "12"="Astrocytes",
                      "13"="MSNs",
                      "14"="Endothelial",
                      "15"="Glutamatergic",
                      "16"="Ependymal",
                      "17"="Immune")

# save the new idents as a column 
brain <- AddMetaData(brain, metadata = Idents(brain), col.name = "broad.labels")

# save the obj ------------------------------------------------------------

saveRDS(brain, file.path(objects_dir,"brain_objects",
                         "05_11_2024_brain_v1.5_broad.labels_added.RDS"))










