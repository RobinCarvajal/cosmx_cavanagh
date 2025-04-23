
####---- SL2 Theme ----####

theme_SL2 <- function () { 
  theme_bw() %+replace% 
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
      plot.background = element_blank(), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      plot.title = element_text(size=12, margin = margin(b = 5),hjust=0,vjust=0.5, face="bold"),
      title = element_text(size = 12, margin = margin(b = 5),hjust=0,vjust=0.5, face="bold"),
      axis.text.y = element_text(size = 10, margin = margin(r = 5),hjust=1,vjust=0.5, face="bold",colour="black"),
      axis.text.x = element_text(size = 10, margin = margin(t = 5),hjust=0.5,vjust=1, face="bold",colour="black"),
      axis.title.y = element_text(size = 11, margin = margin(r = 10),angle = 90,hjust=0.5,vjust=0.5, face="bold"),
      axis.title.x = element_text(size = 11, margin = margin(t = 10),hjust=0.5,vjust=1, face="bold"),
      legend.text=element_text(size=11, face="bold"),
      legend.title=element_blank(), 
      legend.key.size=unit(2.5,"line"),
      plot.margin=unit(c(0.4,0.4,0.4,0.4), "cm")
    )
}


####----- Default GGplot Colours Function -----####

gg_color_hue <- function(n)
{
  # default colour blind palettes
  cblind_palette = c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
  
  # append ggplot colours if groups goes above colour blind palette length
  if (n <= length(cblind_palette))
  {
    return(cblind_palette[1:n])
  }
  else
  {
    nn = n - length(cblind_palette)
    hues = seq(15, 375, length = nn + 1)
    return(append(cblind_palette,hcl(h = hues, l = 65, c = 100)[1:nn]))
  }
}

####---- Save Plot Function -----####

save_plot <- function(ggp,plot_height,plot_width,plot_path)
{
  png(plot_path, height=plot_height, width=plot_width, pointsize=5)
  print(ggp)
  dev.off()
  
  svg(gsub(".png", ".svg", plot_path), height=plot_height/90, width=plot_width/90)
  print(ggp)
  dev.off()
  
  # clears all devices - as a safety measure
  while (dev.cur()>1) dev.off()
}

####----- Clustering Method -----####

cluster_matrix <- function(mx,cluster_x,cluster_y,distance_method,clustering_method,reorder_function)
{
  x.cluster = hclust(Dist(t(mx), method=distance_method), method=clustering_method)
  y.cluster = hclust(Dist(mx, method=distance_method), method=clustering_method)
  
  x.dd = as.dendrogram(x.cluster)
  y.dd = as.dendrogram(y.cluster)
  x.dd.reorder = reorder(x.dd,0,FUN=reorder_function)
  y.dd.reorder = reorder(y.dd,0,FUN=reorder_function)
  x.order = order.dendrogram(x.dd.reorder)
  y.order = order.dendrogram(y.dd.reorder)
  
  if(cluster_x == TRUE && cluster_y == TRUE) 
  {
    mx_clustered = mx[y.order,x.order]
  }
  
  if(cluster_x == TRUE && cluster_y == FALSE)
  {
    mx_clustered = mx[,x.order]
  }
  
  if(cluster_x == FALSE && cluster_y == TRUE)
  {
    mx_clustered = mx[y.order,]
  }
  if(cluster_x == FALSE && cluster_y == FALSE) 
  {
    mx_clustered = mx
  }
  return(mx_clustered)
}

####---- Significant Genes Heatmap Function  ----####

make_significant_genes_heatmap <- function(matrix_scaled,colours,sample_names,cluster_x,cluster_y,distance_method,clustering_method,reorder_function, sample_groupings, default_samples_colours)
{
  
  # tests for enough genes to make a heatmap
  if(nrow(matrix_scaled) >= 2) 
  {
    
    # parse the table
    matrix_scaled_renamed = as.matrix(matrix_scaled)
    colnames(matrix_scaled_renamed) = sample_names
    
    # cluster the table
    matrix_scaled_renamed_clustered = cluster_matrix(matrix_scaled_renamed,cluster_x,cluster_y,distance_method,clustering_method,reorder_function)
    matrix_scaled_renamed_clustered_melted = melt(matrix_scaled_renamed_clustered)
    matrix_scaled_renamed_clustered_melted$X1 = factor(matrix_scaled_renamed_clustered_melted$X1, levels=row.names(matrix_scaled_renamed_clustered))
    matrix_scaled_renamed_clustered_melted$X2 = factor(matrix_scaled_renamed_clustered_melted$X2, levels=colnames(matrix_scaled_renamed_clustered))
    
    # colours
    hm.palette = colorRampPalette(colours)
    
    # make the heatmap
    ggp.hm = ggplot(matrix_scaled_renamed_clustered_melted, aes(x = X2, y = X1, fill = value)) +
      geom_tile() +
      scale_fill_gradientn(colours = hm.palette(100)) +
      ylab('') +
      xlab('') +
      theme_SL2() +
      theme(legend.position="right", legend.title = element_blank(), legend.spacing.x = unit(0.25, 'cm'),axis.text.x = element_text(angle = 90, hjust = 1), axis.ticks=element_blank()) +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0))
    
    
    ## create the rug data
    rug_data = data.frame(default_samples_colours)
    names(rug_data) = "colours"
    row.names(rug_data) = samples
    rug_data$Y = "1"
    rug_data = rug_data[ gsub(" ","_",colnames(matrix_scaled_renamed_clustered)),]
    rug_data$X = row.names(rug_data)
    rug_data$X = factor(row.names(rug_data), row.names(rug_data))
    
    
    # create rug
    ggp.rug = ggplot(rug_data, aes(x = X, y = Y)) +
      geom_tile(fill = rug_data$colours) +
      theme_void() +
      theme(legend.position = "none") +
      ylab('') +
      xlab('') +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_discrete(expand=c(0,0))
    
    # sets the widths
    max_widths = unit.pmax(ggplotGrob(ggp.hm)$widths)
    ggp.rug = ggplotGrob(ggp.rug)
    ggp.rug$widths = max_widths
    
    # combine
    ggp = as.ggplot(grid.arrange(ggp.rug, ggp.hm,  ncol = 1, heights = c(1, 15)))
    
    # clears all devices - as a safety measure
    pdf(NULL)
    while (dev.cur()>1) dev.off()
    
    return(ggp)
  }
  else
  {
    return(ggplot(data.frame()) + theme_SL2() + geom_blank() + ggtitle("There were too few significant genes to plot this."))
  }
}

