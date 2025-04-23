

# FN create Data input from multiple samples -----------------------------------
createData_multisample <- function(obj){
  #parameter# obj: Seurat object (single region, multiple samples)
  #output#    data.input: data input for CellChat
  
  ## Create the metadata from multiple samples 
  # list to store the dataframes
  data_list <- list()
  # loop 
  for (sample in names(obj)){
    data <- obj[[sample]][["RNA"]]$data
    data_list[[sample]] <- data
  }
  
  # number of samples in the object
  n_samples <- length(data_list)
  # standard pipeline
  genes <- lapply(data_list, rownames) # Extracts the rownames of each sample in the data_list
  genes.common <- Reduce(intersect, genes) # Iteratively applies the intersect function across all elements (rownames of each sample) in the list
  
  # Subset and combine data based on genes.common dynamically
  data.input <- do.call(cbind, lapply(data_list, function(sample) sample[genes.common, ]))
  
  return(data.input)
}


# FN Create Metadata from multiple samples -------------------------------------
createMetadata_multisample <- function(obj){
  #parameter# obj: Seurat object (single region, multiple samples)
  
  ## Create the metadata from multiple samples 
  # list to store the dataframes
  meta_list <- list()
  # loop 
  for (sample in names(obj)){
    meta <- obj[[sample]]@meta.data
    meta_list[[sample]] <- data.frame(matrix(ncol = 0, nrow = nrow(meta))) # cc=CellChat
    rownames(meta_list[[sample]]) <- rownames(meta)
    meta_list[[sample]]$samples <- meta$sample_name
    meta_list[[sample]]$broad.labels <- meta$broad.labels
    meta_list[[sample]]$fine.labels <- meta$fine.labels
    # x coords
    meta_list[[sample]]$x_slide_mm <- meta$x_slide_mm
    # y coords
    meta_list[[sample]]$y_slide_mm <- meta$y_slide_mm
  }
  
  # merge the meta_list in a single dataframe
  meta_merged <- do.call(rbind, meta_list)
  # Clean up row names to remove prefixes
  rownames(meta_merged) <- sub("^.*?\\.", "", rownames(meta_merged))
  
  return(meta_merged)
  
}

# FN crate Spatial Locs multisample -------------------------------------------
createSpatialLocs_multisample <- function(obj){
  #parameter# obj: Seurat object (single region, multiple samples)
  
  ## Create the metadata from multiple samples 
  # list to store the dataframes
  spatial_locs_list <- list()
  # loop 
  for (sample in names(obj)){
    coords <- GetTissueCoordinates(obj[[sample]])
    spatial_locs_list[[sample]] <- coords
  }
  
  # merge the meta_list in a single dataframe
  spatial_locs_merged <- do.call(rbind, c(spatial_locs_list, make.row.names=F))
  
  # set cell as rownames
  rownames(spatial_locs_merged) <- spatial_locs_merged$cell
  # delete cell column
  spatial_locs_merged$cell <- NULL
  
  return(spatial_locs_merged)
}

# FN create spatial factors ---------------------------------------------------
createSpatialFactors_multisample <- function(metadata,  conversion.factor = 1000){
  #parameter# meta: metadata with spatial locations
  #parameter# conversion.factor: conversion factor from mm to um (CosMx)
  
  # for multiple samples
  # empty dataframe to store the spatial factors
  spatial.factors <- data.frame()
  # loop
  for (sample in unique(metadata$sample)){
    
    # spatial.locs per sample
    spatial.locs <- metadata %>%
      filter(sample == sample) %>%
      select(x_slide_mm, y_slide_mm)
    
    # standard pipeline
    d <- computeCellDistance(spatial.locs)
    spot.size <- min(d)*conversion.factor # converting the distance
    spatial.factors_sample <- data.frame(ratio = conversion.factor, tol = spot.size/2)
    rownames(spatial.factors_sample) <- sample
    
    # add to the spatial factors dataframe
    spatial.factors <- rbind(spatial.factors, spatial.factors_sample)
    
  }
  
  return(spatial.factors)
}

# FN cellchat drop levels -----------------------------------------------------

cellchat_droplevels <- function(data, meta_cols, clear_idents = TRUE){
  #parameter# data: CellChat object
  #parameter# meta_cols: columns to drop levels
  #parameter# clear_idents: clear idents
  #output#    cleared_data: cleared CellChat object
  
  # dropping levels of idents
  if (clear_idents == TRUE){
    
    data@idents <- droplevels(data@idents)
    
  }
  
  # dropping levels in a loop 
  for (col in meta_cols){
    data@meta[[col]] <- droplevels(data@meta[[col]])
  }
  
  # return the result as a cleared output
  cleared_data <- data
  
  return(cleared_data)
}

# FN standard cellchat pipeline ------------------------------------------------

cellchat_pipeline <- function(data, db, contact_range, interaction_range){
  #parameter# obj: CellChat object
  #parameter# db: CellChat database
  #parameter# contact_range: cell soma/body diameter (mirons)
  #parameter# interaction_range: the maximum interaction/diffusion length of ligands (microns)
  
  # Adding DB
  data@DB <- db
  
  # Subsetting data to include only genes involved in cell interactions
  data <- subsetData(data)
  
  
  # Standard pipeline
  data <- identifyOverExpressedGenes(data)
  print("Over expressed genes identified")
  data <- identifyOverExpressedInteractions(data)
  print("Over expressed interactions identified")

  
  # compute communication probability  
  data <- computeCommunProb(data, type="truncatedMean",trim=0.1,
                            distance.use="FALSE",
                            scale.distance=NULL,
                            contact.dependent=TRUE,
                            interaction.range=interaction_range, 
                            contact.range=contact_range)
  
  print("Communication Probability computed")
  data <- filterCommunication(data, min.cells = 10)
  print("Cummunications filtered")
  data <- computeCommunProbPathway(data)
  print("Communication Probability Pathway computed")
  data <- aggregateNet(data)
  print("Aggregated Cell-Cell communication network calculated")
  
  return(data)
}
