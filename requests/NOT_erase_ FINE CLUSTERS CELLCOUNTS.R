
# immune ------------------------------------------------------------------


for (immune.cluster in c("Immune_0","Immune_1","Immune_2","Immune_3","Immune_4")){
  
  meta <- brain@meta.data
  # subsets per immune cluster
  cluster.meta <- subset(meta, fine.labels == immune.cluster)
  
  table <- table(cluster.meta$sample_name, cluster.meta $condition)
  table_sum <- rbind(table, colSums(table))
  
  print(table_sum)
  
  # save it to a tsv file 
  write.table(table_sum, 
              file = file.path(data_dir,"regional",paste0(immune.cluster,"_sample_condition_table.tsv")), 
              sep = "\t", quote = FALSE)
  
}


# astrocytes --------------------------------------------------------------

astrocyte_names <- c("Astrocytes_0","Astrocytes_1","Astrocytes_2","Astrocytes_3","Astrocytes_4")

for (cluster in astrocyte_names){
  
  meta <- brain@meta.data
  # subsets per cluster
  cluster.meta <- subset(meta, fine.labels == cluster)
  
  table <- table(cluster.meta$sample_name, cluster.meta $condition)
  table_sum <- rbind(table, colSums(table))
  
  print(table_sum)
  
  # save it to a tsv file 
  write.table(table_sum, 
              file = file.path(data_dir,"regional",paste0(cluster,"_sample_condition_table.tsv")), 
              sep = "\t", quote = FALSE)
  
}

