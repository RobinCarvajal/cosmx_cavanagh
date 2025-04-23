# read the object
brain <- readRDS(file.path(data_dir,"objects","brain_objects",
                           "2024_11_05_brain_v1.5_broad.labels_added.RDS"))

library(Nebulosa)
library(viridis)

# substract one sample
sample <- subset(brain, sample_name=="s1056ant_bs3")

# extract spatial coordinates
coords <- sample@meta.data%>%
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
sample@reductions$spatial <- NULL
sample[["spatial"]] <- spatial_red

# checking 
DimPlot(sample, reduction="spatial", group.by = "integrated_clusters", label = TRUE, label.box = TRUE, repel = T)

# plot with nebulosa
genes <- c("Cd3d","Cd3e","Cd3g")
p_list<- Nebulosa::plot_density(sample, reduction = "spatial", 
                                features=genes, 
                                joint = TRUE, combine=FALSE, 
                                size = 0.5)

plot <- p_list[[length(p_list)]]+
  coord_fixed(ratio = 1) 

plot

# save image as an example 
ggsave(file.path(data_dir,"broad_single_cell","plots","nebulosa_combined_spatial_density.png"), plot, width = 10, height = 10, dpi = 300)


# testing normal feature plot 
ImageFeaturePlot(sample, 
                 boundaries="segmentation",
                 features = "Gfap")
