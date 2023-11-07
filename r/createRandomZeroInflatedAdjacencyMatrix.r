createRandomZeroInflatedAdjacencyMatrix <- function(num_vertices, zero_prob) {
  adjacency_matrix <- matrix(0, nrow = num_vertices, ncol = num_vertices)
  
  for (i in 1:(num_vertices - 1)) {
    for (j in (i + 1):num_vertices) {
      if (runif(1) < zero_prob) {
        adjacency_matrix[i, j] <- 0
        adjacency_matrix[j, i] <- 0
      } else {
        weight <- runif(1)
        adjacency_matrix[i, j] <- weight
        adjacency_matrix[j, i] <- weight
      }
    }
  }
  
  return(adjacency_matrix)
}