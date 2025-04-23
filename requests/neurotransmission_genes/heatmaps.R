####---- general libraries ----####

library(ggplot2)
library(ggrepel)
library(reshape)
library(ggridges)

####---- heatmap libraries ----####

library(amap)

####---- table libraries ----####

library(grid)
library(gridExtra)
library(gtable)
library(ggplotify)

####---- DE sample and group names ----####

samples = c("S1055BS1_ANT_C","S1055BS2_ANT_C","S1056BS1_ANT_C","S1056BS2_ANT_C","S1055BS3_ANT_A","S1055BS4_ANT_A","S1056BS3_ANT_A","S1056BS4_ANT_A")
sample_groups = c("ANT_C","ANT_A")
sample_groups = factor(sample_groups, levels = sample_groups)
sample_groupings = c("ANT_C","ANT_C","ANT_C","ANT_C","ANT_A","ANT_A","ANT_A","ANT_A")
sample_groupings = factor(sample_groupings, levels = sample_groups)
samples_by_sample_group = list(c("S1055BS1_ANT_C","S1055BS2_ANT_C","S1056BS1_ANT_C","S1056BS2_ANT_C"),c("S1055BS3_ANT_A","S1055BS4_ANT_A","S1056BS3_ANT_A","S1056BS4_ANT_A"))
comparisons = c("ANT_A vs ANT_C")

####---- sets the working directory. If you move the SL2 folder please update this path ----####

setwd("F:/cosmx_data/regional/differential_expression/condition/anterior/agranular.insular.area/SL2_results/all_genes/de_workflows/ANT_A_vs_ANT_C")

####---- de input files ----####

de_annotated = read.table(file="data/de_annotated.csv", header=TRUE,row.names = 2, sep='\t', quote='',check.names = TRUE)
#de_annotated_sig = read.table(file="data/de_annotated_significant.csv", header=TRUE,row.names = 2, sep='\t', quote='',check.names = TRUE)
neurotransmission_genes=read.csv("F:/cosmx_data/requests/heatmaps/neurotransmission_genes.csv", header=F)
neurotransmission_genes=neurotransmission_genes$V1
de_neurotransmission_genes= de_annotated[neurotransmission_genes,]
####---- de parse ne matrix ----####

ne_matrix = de_annotated[,samples]
#ne_matrix_sig = de_annotated_sig[,samples]
ne_matrix_sig = de_neurotransmission_genes[,samples]

####---- transpose and scale significant genes ne matrix ----####

ne_matrix_sig_scaled = data.frame(t(scale(t(ne_matrix_sig))))
rownames(ne_matrix_sig_scaled) = rownames(ne_matrix_sig)
ne_matrix_sig_scaled[do.call(cbind, lapply(ne_matrix_sig_scaled, is.nan))] <- 0
ne_matrix_sig_scaled = ne_matrix_sig_scaled[is.finite(rowSums(ne_matrix_sig_scaled)), ]


####---- Default Three Tone Heatmap Colours  ----####

default_three_tone_heatmap_colours = c("purple","black","yellow")

####---- Default Sample Labels  ----####

default_sample_labels = gsub("_"," ",c("S1055BS1_ANT_C","S1055BS2_ANT_C","S1056BS1_ANT_C","S1056BS2_ANT_C","S1055BS3_ANT_A","S1055BS4_ANT_A","S1056BS3_ANT_A","S1056BS4_ANT_A")) # note: changing this won't change the order the samples appear on the plots. Merely what they are labelled as.

####---- Default Sample and Sample Group Colours  ----####

number_of_sample_groups = length(sample_groups)
default_sample_group_colours = gg_color_hue(number_of_sample_groups)
default_samples_colours = c(default_sample_group_colours[1],default_sample_group_colours[1],default_sample_group_colours[1],default_sample_group_colours[1],default_sample_group_colours[2],default_sample_group_colours[2],default_sample_group_colours[2],default_sample_group_colours[2])


####---- Significant Genes Heatmap  ----####

plot_height = 2000
plot_width = 500
colours = default_three_tone_heatmap_colours
sample_names = default_sample_labels
distance_method = "spearman"
clustering_method = "average"
reorder_function = "average"

cluster_x = FALSE
cluster_y = TRUE

out = paste0("F:/cosmx_data/requests/heatmaps/significant_genes_heatmap_", "anterior_cingulate.png")

ggp = make_significant_genes_heatmap(ne_matrix_sig_scaled,colours,sample_names,cluster_x,cluster_y,distance_method,clustering_method,reorder_function, sample_groupings, default_samples_colours)
save_plot(ggp,plot_height,plot_width,out)

cluster_x = TRUE
cluster_y = TRUE

# ggp = make_significant_genes_heatmap(ne_matrix_sig_scaled,colours,sample_names,cluster_x,cluster_y,distance_method,clustering_method,reorder_function, sample_groupings, default_samples_colours)
# save_plot(ggp,plot_height,plot_width,"plots/significant_genes_heatmap/significant_genes_heatmap_column_clustered.png")

