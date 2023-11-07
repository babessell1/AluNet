createAdjacencyMatrixFromTSV <- function(file_path, chromosome) {
  # Initialize variables
  adjacency_matrix <- NULL
  vertices <- character(0)
  num_vertices <- 0

  # Read the TSV file line by line
  con <- file(file_path, "r")
  while (length(line <- readLines(con, n = 1)) > 0) {
    fields <- strsplit(line, "\t")[[1]]

    # Check if the line belongs to the specified chromosome
    fragment1 <- fields[1]
    fragment2 <- fields[2]

    if (startsWith(fragment1, chromosome) && startsWith(fragment2, chromosome)) {
      if (!(fragment1 %in% vertices)) {
        vertices <- c(vertices, fragment1)
        num_vertices <- num_vertices + 1
      }
      if (!(fragment2 %in% vertices)) {
        vertices <- c(vertices, fragment2)
        num_vertices <- num_vertices + 1
      }
    }
  }
  close(con)

  # Initialize the adjacency matrix based on the number of vertices
  adjacency_matrix <- matrix(0, nrow = num_vertices, ncol = num_vertices)
  rownames(adjacency_matrix) <- colnames(adjacency_matrix) <- vertices

  # Read the file again and assign weights
  con <- file(file_path, "r")
  while (length(line <- readLines(con, n = 1)) > 0) {
    fields <- strsplit(line, "\t")[[1]]

    # Check if the line belongs to the specified chromosome
    fragment1 <- fields[1]
    fragment2 <- fields[2]
    raw_count <- as.numeric(fields[4])

    if (startsWith(fragment1, chromosome) && startsWith(fragment2, chromosome)) {
      adjacency_matrix[fragment1, fragment2] <- raw_count
      adjacency_matrix[fragment2, fragment1] <- raw_count
    }
  }
  close(con)

  return(adjacency_matrix)
}
