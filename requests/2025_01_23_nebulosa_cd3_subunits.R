
# load brain obj 
brain <- readRDS(file.path(objects_dir,"brain_objects",
                         "2025_01_22_brain_v1.7_region.labels_added.RDS"))

# prepare for nebulosa plots

prepareNebulosa <- function(obj)
{
  # extract spatial coordinates
  coords <- obj@meta.data%>%
    dplyr::select(x_slide_mm,y_slide_mm)%>%
    as.matrix()
  # renaming cols
  colnames(coords) <- c("SPATIAL_1", "SPATIAL_2")
  # creating a reduction object
  spatial_red <- CreateDimReducObject(
    embeddings = coords, # Matrix for cell embeddings
    key = "SPATIAL_",               # Prefix for dimension names
    assay = "RNA"              # Assay the dim reduction applies to
  )
  # add the UMAP coords to the seurat object
  obj@reductions$spatial <- NULL
  obj[["spatial"]] <- spatial_red
  
  return(obj)
  
}

# plot nebulosa spatial fucntion --------------------------------------
DimPlot(brain, reduction="spatial", group.by = "integrated_clusters", label = TRUE, label.box = TRUE, repel = T)

plotNebulosa_spatial <- function(obj, sample, genes)
{
  
  sample.obj <- subset(obj, sample_name == sample)
  
  # plot with nebulosa
  p_list<- Nebulosa::plot_density(sample.obj, reduction = "spatial", 
                                  features=genes, 
                                  joint = TRUE, combine=FALSE, 
                                  size = 0.5)
  
  plot <- p_list[[length(p_list)]]+
    coord_fixed(ratio = 1) 
  
  return(plot)
  
}

# generating the plots --------------------------------------
brain <- prepareNebulosa(brain)

genes <- c("Cd3d","Cd3e","Cd3g")
genes <- "Cd3d"

nebulosa_plots <- list()

for (sample in unique(brain$sample_name))
{
  plot <- plotNebulosa_spatial(brain, sample = sample, genes = genes)
  nebulosa_plots[[sample]] <- plot
}

# save plots --------------------------------------

# create a folder to store the plots
plots_dir <- file.path(data_dir, "requests", "nebulosa_cd3_subunits")
dir.create(plots_dir)

# save in a loop 
for (sample in unique(brain$sample_name))
{
  plot <- nebulosa_plots[[sample]]
  ggsave(file.path(plots_dir, paste0(sample,"_cd3_subunits.png")), plot, width = 10, height = 10)
}


# put two together, same scale --------------------------------------

# get the list of plots
