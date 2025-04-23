# run needed libraries
source("F:/cosmx_scripts/Set_up_HOME.R")

# Load the processed cellchat objects ------------------------------------------
cellchat_objects_processed <- readRDS(file = file.path(data_dir, "regional","cellchat_objects_processed.RDS"))

# Extracting interaction tables ------------------------------------------------
cellchat_interaction_tables <- list()
# loop
for (region in names(cellchat_objects_processed)){
  data <- cellchat_objects_processed[[region]]
  interactions <- subsetCommunication(data)
  cellchat_interaction_tables[[region]] <- interactions
}


data <- cellchat_objects_processed[["caudoputamen"]]

groupSize <- as.numeric(table(data@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(data@net$count, vertex.weight = rowSums(data@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(data@net$weight, vertex.weight = rowSums(data@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

pathways.show<-c("CSF","PTN")
pathways.show <- (unique(cellchat_interaction_tables[["anterior_cingulate_cortex"]]$pathway_name))
netAnalysis_contribution(data, signaling = pathways.show)

spatialFeaturePlot(data, sample.use="s1055ant_bs4", features="Psap")

netVisual_heatmap(data, measure = "count", color.heatmap = "Blues")

netVisual_heatmap(data, measure = "weight", color.heatmap = "Blues")

# visualisation of cell-celll communication netweork 
pathways.show <- c("SOMATOSTATIN") 
# Circle plot
par(mfrow=c(1,1), xpd=TRUE)
netVisual_aggregate(data, signaling = pathways.show, layout = "circle")

# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(data, 
                    signaling = pathways.show, 
                    sample.use = "s1055ant_bs4", layout = "spatial", 
                    edge.width.max = 2, vertex.size.max = 1,
                    alpha.image = 0.2, vertex.label.cex = 0)

# Compute the network centrality scores
data <- netAnalysis_computeCentrality(data, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(data, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# USER can show this information on the spatial transcriptomics when visualizing a signaling network, e.g., bigger circle indicates larger incoming signaling
par(mfrow=c(1,1))
netVisual_aggregate(data, signaling = pathways.show, sample.use = "s1055ant_bs1", 
                    layout = "spatial", edge.width.max = 2, alpha.image = 0.2, 
                    vertex.weight = "incoming", vertex.size.max = 6, vertex.label.cex = 0)

# compute controbutiuon of each L-R pair to the overall signalling pathway 
netAnalysis_contribution(data, signaling = pathways.show)

# Take an input of a few genes
spatialFeaturePlot(data, features = c("Sst","Sstr4"), sample.use = "s1055ant_bs4", 
                   point.size = 0.8, color.heatmap = "Reds", direction = 1)

# Take an input of a ligand-receptor pair
spatialFeaturePlot(data, pairLR.use = "SST_SSTR4", sample.use = "s1055ant_bs4", 
                   point.size = 0.5, do.binary = FALSE, cutoff = 0.05, 
                   enriched.only = F, color.heatmap = "Reds", direction = 1)

# Take an input of a ligand-receptor pair and show expression in binary
spatialFeaturePlot(data, pairLR.use = "SST_SSTR4", sample.use = "s1055ant_bs4",
                   point.size = 1.5, do.binary = TRUE, cutoff = 0.05, 
                   enriched.only = F, color.heatmap = "Reds", direction = 1)




##########################
# calcualting contact range
locs <- data@images$coordinates
d <- computeCellDistance(locs)
contact.range <- mean(d)
