#' Plot the output of the Leiden algorithm, coloring nodes based on community
#' 
#' @param result The result of the runLeiden function.
#' @return A plot of the graph with nodes colored based on community.
#' @export
plotLeiden <- function(result, list, Alus = FALSE) {
  gdf <- graph.data.frame(result$graph)
  
  # Get node indices from the graph
  node_indices = V(gdf)$name
  
  # Create map of node names to their communities
  node_community_map <- setNames(result$communities, result$node_names)
  
  # Assign membership in graph dataframe based on result of Leiden
  V(gdf)$membership <- node_community_map[node_indices]
  
  # Create colormap for the communities
  num_communities <- length(unique(result$communities))
  community_colors <- rainbow(num_communities)
  colors <- community_colors[as.factor(V(gdf)$membership)]
  
  # Define the edge colors based on ALU group
  # alu_colors <- c("AluJ" = "purple", "AluS" = "black", "AluY" = "red")
  # edge_alpha <- 0.8
  
  # Matching edges from the graph to the list and assigning colors
  # matched_indices <- match(paste(E(gdf)$from, E(gdf)$to), paste(list$from, list$to))
  # print(matched_indices)
  # default_color <- "grey"  # Default color for unmatched edges
  # edge_colors <- rep(default_color, ecount(gdf))
  # valid_indices <- !is.na(matched_indices)
  # edge_colors[valid_indices] <- adjustcolor(alu_colors[list$group[matched_indices[valid_indices]]], alpha.f = edge_alpha)
  # E(gdf)$color <- edge_colors
  
  # Plot logic
  plot_logic <- function() {
    plot(
      gdf,
      layout = layout_with_fr(gdf),
      vertex.color = colors,
      vertex.size = 4,
      vertex.frame.color = colors,
      edge.color = E(gdf)$color,  # Edge colors
      edge.width = E(gdf)$weight * 8,  # Edge width
      vertex.label = ifelse(Alus, NA, V(gdf)$name),
      edge.arrow.size = 0,
      main = ifelse(Alus, "Candidate Chromatin Modulating Alu Communities (Chr22)", "Our Rcpp-based Clustering")
    )
    if (!Alus) {
      legend("topright", legend = unique(V(gdf)$membership),
             col = community_colors, pch = 19,
             title = "Communities", cex = 0.8)
    }
  }
  
  # Execute the plot logic
  plot_logic()
}
