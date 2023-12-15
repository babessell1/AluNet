

#' Read DFAM Data and Filter ALU Sequences
#'
#' Reads a DFAM data file in chunks, filters for sequences starting with 'ALU', and saves them to CSV files.
#'
#' @param file_path Path to the DFAM data file.
#' @param output_dir Directory where the filtered CSV files will be saved.
#' @return NULL
#' @import data.table
#' @importFrom stringr grepl
read_dfam <- function(file_path, output_dir){
  # Set parameters
  output_dir <- paste(output_dir, "selection_by_chr", sep = "/")
  # Create the directory
  if (!file.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  chunksize <- 10^6
  skip_rows <- 0 
  
  # Function to check if a string starts with 'ALU'
  startsWithAlu <- function(x) {
    grepl("^ALU", toupper(x))
  }
  
  header <- colnames(fread(file_path, nrows = 1))
  print(header)
  # Read and process the file in chunks
  continue <- TRUE
  while(continue) {
    tryCatch({
      # Read a chunk of the file
      dt <- fread(file_path, skip = skip_rows, nrows = chunksize, header = FALSE, col.names = header)
     
      if (nrow(dt) == 0) {
        print(paste("Reached end of file or no more data from skip_rows =", skip_rows))
        continue = FALSE
        return()
      }
      
      # Filter rows where family_name starts with 'ALU'
      alu_rows <- dt[startsWithAlu(dt$family_name), ]
    
      names(alu_rows)[names(alu_rows) == "#seq_name"] <- "seq_name"
    
    # Write to CSV
    fwrite(alu_rows, paste0(output_dir, "/alu_data", skip_rows, ".csv"), col.names = TRUE)
    print(paste("We have finished the", skip_rows, "to", (skip_rows + chunksize)))
    
    
    # Update skip_rows for the next iteration
    skip_rows <- skip_rows + chunksize
    }, error = function(e) {
      # Handle the error
      print(paste("Something wrong happens", skip_rows, ": ", e$message))
      if(grepl("skip=", e$message)) {
        print("Reached past the end of the file. Stopping.")
        continue <<- FALSE
      }
    })
  }
}

#' Categorize ALU Sequences by Chromosome
#'
#' Reads ALU sequence data and categorizes them by chromosomes based on a given Hi-C data file.
#'
#' @param hic_file_path Path to the Hi-C data file.
#' @param output_dir Directory where the categorized CSV files will be saved.
#' @return Invisible NULL. The function is used for its side effects of reading, filtering, and writing data.
#' @import dplyr
#' @import purrr
#' @import readr
categorize_dfam <- function(hic_file_path, output_dir){
  running_dir <- paste(output_dir, "selection_by_chr", sep = "/")
  # place all your data to the uncleaned data folder
  output_path_cleaned <- paste0(output_dir, "/selection_by_chr_cleaned")
  # Create the directory
  if (!file.exists(output_path_cleaned)) {
    dir.create(output_path_cleaned)
  }
  
  hic_data <- read.delim(hic_file_path, header = TRUE)
  hic_data$chromosome <- sub("^([^.]*)\\..*", "\\1", hic_data$frag1)
  unique_chromosomes <- unique(hic_data$chromosome)
  
  read_and_filter <- function(file, chromosome) {
    data <- read.csv(file)
    filtered_data <- data %>% filter(seq_name == chromosome)
    return(filtered_data)
  }
  
  files <- list.files(running_dir, pattern="*.csv", full.names = TRUE)
  # we first extract the chromosome's alu elements
  # print(files)
  for(i in unique_chromosomes){
    cat(paste0("now we are working on chromosome: ", i))
    
    merged_data <- map_df(files,  ~read_and_filter(., i))
    
    print(nrow(merged_data))
    
    output_file <- paste0(output_path_cleaned, "/", i, ".csv")
    write.csv(merged_data, output_file)
  }
}

#' Select ALU Sequences with the Highest Probability
#'
#' For each position in the data frame, selects the ALU sequences with the highest probability.
#'
#' @param output_dir Directory containing cleaned ALU sequence data.
#' @return Invisible NULL. The function is used for its side effects of selecting and writing data.
#' @import data.table
select_the_highest_probability_alu <- function(output_dir){
  # we first extract the key alu elements, which with highest probability
  # for each position in the data frame
  # then we merge the data frame with the whole data.
  current_directory <- paste0(output_dir, "/selection_by_chr_cleaned")
  output_path <- paste0(output_dir, "/highest_probability")
  # Create the directory
  if (!file.exists(output_path)) {
    dir.create(output_path)
  }
  
  files <- list.files(current_directory, pattern="*", full.names = TRUE)
  
  for(i in files){
    last_item <- basename(i)
    
    cat(paste0("now we are working on chromosome: ", last_item))
    
    test_file <- read.csv(i)
    test_file <- data.table(test_file)
    
    df_selected <- test_file[, .SD[which.min(e.value)], by = .(ali.st, ali.en)]
    print(nrow(df_selected))
    write.csv(df_selected, paste0(output_path, '/', last_item))
  }
}

#' Main Function to Process ALU Data
#'
#' This is a pipeline function that sequentially processes ALU data from a DFAM file. 
#' It first reads and filters ALU sequences from a given DFAM file (read_dfam), 
#' then categorizes these sequences by chromosomes based on a Hi-C data file (categorize_dfam), 
#' and finally selects the ALU sequences with the highest probability for each position in the data frame 
#' (select_the_highest_probability_alu).
#'
#' @param hic_file_path Path to the Hi-C data file, used in categorizing ALU sequences by chromosome.
#' @param file_path Path to the DFAM data file, which contains ALU sequences to be processed.
#' @param output_dir Base directory where the processed output data will be saved.
#' 
#' @details 
#' The `read_dfam` function reads a DFAM file in chunks, filters for sequences starting with 'ALU', and saves them to CSV files.
#' The `categorize_dfam` function reads the filtered ALU data and categorizes them by chromosomes based on the Hi-C data file.
#' The `select_the_highest_probability_alu` function selects the ALU sequences with the highest probability for each position.
#'
#' @return path of the folder where the final output is stored.
clean_alu <- function(hic_file_path, file_path, output_dir) {
  read_dfam(file_path, output_dir)
  categorize_dfam(hic_file_path, output_dir)
  select_the_highest_probability_alu(output_dir)
  return (paste0(output_dir, "/highest_probability"))
}