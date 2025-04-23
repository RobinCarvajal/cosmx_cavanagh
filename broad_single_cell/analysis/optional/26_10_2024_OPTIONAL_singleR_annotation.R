
# load object 
brain <- readRDS(file.path(main_dir, "New_Analysis","objects",
                           "26_10_2024_v1.0_brain_integrated.RDS"))

# call singleR library 
library(SingleR)
library(celldex) # this is used for the reference
library(scRNAseq) # this is used for the scRNAseq data
#library(pheatmap) # this is used for the heatmap


# singleR manual  ---------------------------------------------------------

# https://bioconductor.org/books/devel/SingleRBook/sc-mode.html

# reference ---------------------------------------------------------------

# bulk RNAseq ref

general_ref <- celldex::MouseRNAseqData() # using this one 
# immune_ref <- celldex::ImmGenData() 

# single cell ref # TEST SOME DAY

# (manual: https://bioconductor.org/packages/release/data/experiment/manuals/scRNAseq/man/scRNAseq.pdf)
# campbell_single_cell_ref <- scRNAseq::CampbellBrainData() #2017
# chen_single_cell_ref <- scRNAseq::ChenBrainData() #2017
# jessa_single_cell_ref <- scRNAseq::JessaBrainData() #2019
# marques_single_cell_ref <- scRNAseq::MarquesBrainData() #2016
# romanov_single_cell_ref <- scRNAseq::RomanovBrainData() #2017
# tasic_single_cell_ref <- scRNAseq::TasicBrainData() #2015
# usoskin_ref <- scRNAseq::UsoskinBrainData() #2015 CRAP
# lamanoo_ref <- scRNAseq::LaMannoBrainData("mouse-adult") #2016 promising


# general ref (fine) -------------------------------------------------------------

# Extracting brain counts
counts <- LayerData(brain, assay = "RNA", layer = "counts")

# Predictions
predictions.fine <- SingleR(test = counts, 
                       ref = general_ref, 
                       labels = general_ref$label.fine # using fine/detailed labels
                       )

# Get the frequency table
predictions_table.fine <- as.data.frame(table(predictions.fine$labels))

# assign labels to sample_obj
brain$singleR_labels.fine <- predictions.main$labels[match(rownames(brain@meta.data), 
                                                           rownames(predictions.fine))]

# visualise the new labels
DimPlot(brain, 
        reduction = 'umap.integrated', 
        group.by = c("integrated_clusters",'singleR_labels.fine'), 
        label = TRUE)

# general ref (main) ------------------------------------------------------

# predictions
predictions.main <- SingleR(test = counts, 
                       ref = general_ref, 
                       labels = general_ref$label.main # using fine/detailed labels
                       )

# Get the frequency table
predictions_table.main <- as.data.frame(table(predictions.main$labels))

# assign labels to sample_obj
brain$singleR_labels.main <- predictions.main$labels[match(rownames(brain@meta.data), 
                                                           rownames(predictions.main))]

# visualise the new labels
DimPlot(brain, 
        reduction = 'umap.integrated', 
        group.by = c("integrated_clusters",'singleR_labels.main'), 
        label = TRUE)


# save the object ---------------------------------------------------------

saveRDS(brain, file.path(main_dir, "New_Analysis","objects","26_10_2024_v1.0_brain_SingleR.RDS"))


# Optional Diagnostics ----------------------------------------------------

# # Score Heatmap
# score_heatmap_pl <- plotScoreHeatmap(predictions)
# 
# save_plot(score_heatmap_pl, 
#           file = file.path(sample_output_SingleR, "score_heatmap_pl.png"), 
#           w=2000, h=1500)
# 
# # Delta Distribution
# delta_dist_pl <- plotDeltaDistribution(predictions)
# 
# save_plot(delta_dist_pl, 
#           file = file.path(sample_output_SingleR, "delta_dist_pl.png"), 
#           w=2000, h=1500)
# 
# 
# # Score Distribution
# tab <- table(Assigned=predictions$labels, 
#              Clusters=sample_obj$harmony_clusters)
# 
# pheatmap_pl <- pheatmap(log10(tab+10), color = colorRampPalette(c('white','blue'))(10))





