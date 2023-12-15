#' Plot the output of the Leiden algorithm, coloring nodes based on community
#' 
#' @param result The result of the runLeiden function.
#' @return A plot of the graph with nodes colored based on community.
#' @export
plotLeiden <- function(result) {
    gdf <- graph.data.frame(result$graph)

    # get node indices from the graph
    node_indices = V(gdf)$name

    # create map of node names to their communitys
    node_community_map <- setNames(result$communities, result$node_names)

    # assign membership in graph dataframe based on result of Leiden
    V(gdf)$membership <- node_community_map[node_indices]

    # create colormap for the communities
    get_community_colors <- function(num_communities) {
    rainbow(num_communities)
    }

    # map communities to the colormap
    num_communities <- length(unique(result$communities))
    community_colors <- get_community_colors(num_communities)
    colors <- community_colors[as.factor(V(gdf)$membership)]

    # Plot the graph
    plot(
    gdf,
    layout = layout_with_fr(gdf),
    vertex.color = colors,  # Color nodes based on community
    vertex.frame.color = colors,
    #edge.color = color_palette[cut(E(gdf)$weight, breaks = 100)],  # Map edge color from cold to hot
    vertex.label = V(gdf)$name,
    edge.arrow.size = 0,  # Set arrow size to 0 to remove arrows
    main = "Our Rcpp-based Clustering"
    )

    # Add a legend for community colors
    legend("topright", legend = unique(V(gdf)$membership),
        col = community_colors, pch = 19,
        title = "Communities", cex = 0.8)

}

