graphFromFitHIC <- function(tsv_file, chromosome) {
  # Read the TSV file into a data frame
  data <- read.table(tsv_file, header = TRUE, sep = "\t")

  # Filter data for the specified chromosome
  filtered_data <- subset(data, grepl(chromosome, frag1) | grepl(chromosome, frag2))

  # Create a list to store edge information
  edges_list <- list(
    from = as.vector(filtered_data$frag1),
    to = as.vector(filtered_data$frag2),
    weight = as.numeric(as.vector(filtered_data$rawCount))
  )

  return(edges_list)
}