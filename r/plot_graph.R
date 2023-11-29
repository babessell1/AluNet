plot.graph <- function(g_list, filename, small_nodes=FALSE) {

    graph_data <- createGraphFromList(g_list)
    edges <- g_list$edges
    weights <- g_list$weights
    #print(graph_data)
    g <- graph.data.frame(graph_data)
    #print(g)
    edge_weights <- E(g)$w
    V(g)$label <- NA
    E(g)$label <- NA
    layout <- layout_with_fr(g, dim=2)

    # Increase the size of the plot
    par(mar = c(1, 1, 1, 1))  # Adjust margins if needed
    layout <- layout_with_fr(g, dim = 2)

    # Set the size of the plot device using the 'width' and 'height' parameters
    png(filename, width = 2000, height = 2000)  # Adjust width and height as needed

    if (small_nodes) {    
        # Plot the graph with smaller nodes
        plot.igraph(g, layout = layout, vertex.size = .01, vertex.label.dist = 2, edge.alpha = 1, edge.arrow.size = 0)

        dev.off()  # Close the plot device
        plot.igraph(g, layout = layout, vertex.size = .01, vertex.label.dist = 2, edge.alpha = 1, edge.arrow.size = 0)
    } else {
        # Plot the graph with larger nodes
        plot.igraph(g, layout = layout, vertex.size = 1, vertex.label.dist = 2, edge.alpha = 1, edge.arrow.size = 0)

        dev.off()  # Close the plot device
        plot.igraph(g, layout = layout, vertex.size = 1, vertex.label.dist = 2, edge.alpha = 1, edge.arrow.size = 0)
    }
}