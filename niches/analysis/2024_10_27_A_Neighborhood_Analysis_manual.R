# libraries
library(Seurat)
#library(here)
library(ggplot2)
library(dplyr)

# load brain object 
brain <- readRDS(file.path(main_dir, "New_Analysis","objects",
                           "26_10_2024_brain_integrated.RDS"))

# get data from one sample
s1055ant_bs1 <- subset(brain, subset = sample_name == "s1055ant_bs1")
meta <- s1055ant_bs1@meta.data

#get the centroids matrix
matrix <- meta[,c("x_slide_mm","y_slide_mm")]

library(dbscan)

radious <- 0.05 # 0.05 mm is 50 micrometers

nn <- frNN(x=matrix, eps=radious)

matrix[nn$id$s1055ant_bs1_c_1_12_390,] #random cell


# neighborhood analysis ---------------------------------------------------


all.equal(colnames(s1055ant_bs1), names(nn$id))


# this takes me 7mins, can figure a better way to do it
x <- purrr::map(nn$id, ~s1055ant_bs1$integrated_clusters[.x] %>% table())

nn_matrix <- do.call(rbind,x)

nrow(nn_matrix) #36357 rows

# create a seurat object 

nn_obj <- CreateSeuratObject(counts = t(nn_matrix),  min.features = 5) #  try with 0 later

nn_obj <- SCTransform(nn_obj, vst.flavor = "v2")

nn_obj <- RunPCA(nn_obj, npcs = 30, features = rownames(nn_obj))
ElbowPlot(nn_obj)

nn_obj <- FindNeighbors(nn_obj, reduction = "pca", dims = 1:10)
nn_obj <- FindClusters(nn_obj, resolution = 0.7)

nn_obj <- RunUMAP(nn_obj, dims = 1:9)
DimPlot(nn_obj)

# put the neighborhood back into the object

old_meta<- s1055ant_bs1@meta.data %>% 
  tibble::rownames_to_column(var= "cell_ident")

nn_meta<- nn_obj@meta.data %>%
  tibble::rownames_to_column(var= "cell_ident") %>%
  select(cell_ident, SCT_snn_res.0.7)

## note, we filtered out some cells for the neighborhood analysis
new_meta<- old_meta %>%
  left_join(nn_meta)

new_meta<- as.data.frame(new_meta)
rownames(new_meta)<- old_meta$cell_ident

s1055ant_bs1@meta.data<- new_meta

# replace NA in SCT_snn_res.0.7 with "FAIL"
s1055ant_bs1$SCT_snn_res.0.7 <- as.character(s1055ant_bs1$SCT_snn_res.0.7)
s1055ant_bs1$SCT_snn_res.0.7[is.na(s1055ant_bs1$SCT_snn_res.0.7)] <- "FAIL"
s1055ant_bs1$SCT_snn_res.0.7 <- as.factor(s1055ant_bs1$SCT_snn_res.0.7)
levels(s1055ant_bs1$SCT_snn_res.0.7)

# remove values that are FAIL in SCT_snn_res.0.7
temp_sample <- subset(s1055ant_bs1, subset = SCT_snn_res.0.7 != "FAIL")

# change default to segmentation 
DefaultBoundary(temp_sample[["s1055ant"]]) <- "segmentation"

## the cells are colored by the clustering of the cells by expression
p1<- ImageDimPlot(temp_sample, fov = "s1055ant", cols = "polychrome", axes = TRUE,
                  border.size = 0.01)


## the cells are colored by the clustering of the cells by neighborhood 
p2<- ImageDimPlot(temp_sample, fov = "s1055ant", cols = "polychrome", axes = TRUE, 
                  group.by = "SCT_snn_res.0.7",
                  border.size = 0.01)

p1 + p2



# save p2 plot
ggsave(file.path(main_dir, "New_Analysis","plots",
                 "27_10_2024_s1055ant_bs1_neighborhood.png"), p2, width = 10, height = 10, dpi = 600)

# save p1
ggsave(file.path(main_dir, "New_Analysis","plots",
                 "27_10_2024_s1055ant_bs1_expression.png"), p1, width = 10, height = 10, dpi = 600)


# which integrated cclusters are within SCT_snn_res.0.7
table(s1055ant_bs1$integrated_clusters, s1055ant_bs1$SCT_snn_res.0.7)

