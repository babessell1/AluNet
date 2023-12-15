#' Clean ALU Data from DFAM
#'
#' This function processes ALU data from the DFAM database. It reads the input file in chunks,
#' filters ALU sequences, writes them to CSV files, categorizes these sequences by chromosome,
#' and selects the ALU elements with the highest probability (minimum e-value) for each position.
#' The function creates several directories for intermediate and final outputs.
#'
#' @param input_file_path The path to the input file containing ALU data.
#' @param output_directory The base directory where output files and folders will be created.
#' @param chunksize The number of rows to read at a time; defaults to 1,000,000.
#' @import data.table
#' @import dplyr
#' @export
clean_alu <- function(input_file_path, output_directory, chunksize = 10^6) {
  library(data.table)
  library(dplyr)
  
  # Check if output directories exist and create them if not
  output_dir_selection_by_chr <- file.path(output_directory, "selection by chr")
  output_dir_selection_by_chr_cleaned <- file.path(output_directory, "selection by chr cleaned")
  output_dir_highest_probability <- file.path(output_directory, "highest probability")
  
  dirs_to_create <- c(output_dir_selection_by_chr, output_dir_selection_by_chr_cleaned, output_dir_highest_probability)
  
  for(dir in dirs_to_create) {
    if (!dir.exists(dir)) {
      dir.create(dir)
    }
  }
  
  # Function to check if a string starts with 'ALU'
  startsWithAlu <- function(x) {
    grepl("^ALU", toupper(x))
  }
  
  # Read and process the input file in chunks
  skip_rows <- 0
  header <- colnames(fread(input_file_path, nrows = 1))
  
  while(TRUE) {
    dt <- fread(input_file_path, skip = skip_rows, nrows = chunksize, header = FALSE, col.names = header)
    if (nrow(dt) == 0) break
    
    alu_rows <- dt[startsWithAlu(dt$family_name), ]
    names(alu_rows)[names(alu_rows) == "#seq_name"] <- "seq_name"
    fwrite(alu_rows, file.path(output_dir_selection_by_chr, paste0("alu_data", skip_rows, ".csv")), col.names = TRUE)
    cat(paste("We have finished the", skip_rows, "to", (skip_rows + chunksize)))
    skip_rows <- skip_rows + chunksize
  }
  
  # Process each chromosome
  hic_data <- fread(file.path(output_directory, "hic_data.txt"))
  hic_data$chromosome <- sub("^([^.]*)\\..*", "\\1", hic_data$frag1)
  unique_chromosomes <- unique(hic_data$chromosome)
  
  read_and_filter <- function(file, chromosome) {
    data <- fread(file)
    data[seq_name == chromosome, ]
  }
  
  files <- list.files(output_dir_selection_by_chr, pattern="*.csv", full.names = TRUE)
  
  for(chromo in unique_chromosomes) {
    merged_data <- purrr::map_df(files, ~read_and_filter(., chromo))
    fwrite(merged_data, file.path(output_dir_selection_by_chr_cleaned, paste0(chromo, ".csv")))
  }
  
  # Process for highest probability
  files <- list.files(output_dir_selection_by_chr_cleaned, pattern="*.csv", full.names = TRUE)
  
  for(file_path in files) {
    chromosome_name <- basename(file_path)
    dt <- fread(file_path)
    
    cat(paste0("now we are working on chromosome: ", chromosome_name))
    
    selected_data <- dt[, .SD[which.min(e.value)], by = .(ali.st, ali.en)]
    fwrite(selected_data, file.path(output_dir_highest_probability, chromosome_name))
  }
}
