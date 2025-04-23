
# load aldara cellchat object
aldara_cellchat_objects_processed <- readRDS(file.path(data_dir,"regional",
                                                       "aldara_cellchat_objects_processed.RDS"))

# load aldara cellchat object
control_cellchat_objects_processed <- readRDS(file.path(data_dir,"regional",
                                                        "control_cellchat_objects_processed.RDS"))

# get thalamus for aldara and control
aldara_thalamus <- aldara_cellchat_objects_processed[["thalamus"]]
control_thalamus <- control_cellchat_objects_processed[["thalamus"]]

# merge them in a single object
thalamus_list <- list(C=control_thalamus, A=aldara_thalamus)
object.list <- thalamus_list
thalamus_merged <- mergeCellChat(thalamus_list, add.names = names(thalamus_list))
cellchat <- thalamus_merged


# plot1 
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg1

# plot2 
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# plot3
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

# plot4
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# plot5
library(plyr)
group.cellType_old <- levels(control_thalamus@idents)
group.cellType_new <- c("debris","debris","neuron","neuron","vascular","glia","glia","neuron","immune","vascular","vascular")
group.cellType <- mapvalues(group.cellType_old, from = group.cellType_old, to = group.cellType_new)
group.cellType <- factor(group.cellType, levels = unique(group.cellType))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})


weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

# plot6
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

# plot7

for (object in names(object.list)){
  object.list[[object]] <- netAnalysis_computeCentrality(object.list[[object]])
}
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


# plot8

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Astrocytes")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Immune")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2))


# plot9

# this needs this 
reticulate::py_install(packages = 'umap-learn')

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

# plot10

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)


# plot11
rankSimilarity(cellchat, type = "functional")
#> Compute the distance of signaling networks between datasets 1 2

# left here





