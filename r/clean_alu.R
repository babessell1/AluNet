#' Clean ALU Data from DFAM
#'
#' This function performs several steps to clean ALU data from DFAM. It reads the 'hg38.hits' file in chunks,
#' filters ALU sequences, writes them to CSV files, categorizes these sequences by chromosome,
#' and finally selects the ALU elements with the highest probability (minimum e-value) for each position.
#' NOTE: this function should run pretty long time (hours)
#'
#' The function assumes a specific directory structure and creates necessary directories if they do not exist.
#'
#' @return NULL. The function operates through side effects: reading, filtering,
#' categorizing, selecting, and writing data to files. It does not return any value.
#'
#' @examples
#' clean_alu()
#'
#' @export
clean_alu <- function(){
  file_path <- getwd()
  download_directory <- paste(file_path, "uncleaned data", sep = "/")
  file_name <- "hg38.hits"
  hg38_file <- paste(download_directory, file_name, sep = "/")
  library(data.table)
  # Set parameters
  output_dir <- paste(file_path, "selection by chr", sep = "/")
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
  
  header <- colnames(fread(hg38_file, nrows = 1))
  print(header)
  # Read and process the file in chunks
  while(TRUE) {
    # Read a chunk of the file
    dt <- fread(hg38_file, skip = skip_rows, nrows = chunksize, header = FALSE, col.names = header)
    # If the chunk is empty, break the loop
    if (nrow(dt) == 0) {
      break
    }
    
    # Filter rows where family_name starts with 'ALU'
    alu_rows <- dt[startsWithAlu(dt$family_name), ]
    
    names(alu_rows)[names(alu_rows) == "#seq_name"] <- "seq_name"
    
    # Write to CSV
    fwrite(alu_rows, paste0(output_dir, "/alu_data", skip_rows, ".csv"), col.names = TRUE)
    print(paste("We have finished the", skip_rows, "to", (skip_rows + chunksize)))
    
    
    # Update skip_rows for the next iteration
    skip_rows <- skip_rows + chunksize
  }
  
  running_dir <- paste(file_path, "selection by chr", sep = "/")
  # place all your data to the uncleaned data folder
  output_path_cleaned <- paste0(file_path, "/selection by chr cleaned")
  # Create the directory
  if (!file.exists(output_path_cleaned)) {
    dir.create(output_path_cleaned)
  }
  
  hic_data <- read.delim(paste0(file_path, "/hic_data.txt"), header = TRUE)
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
  
  # we first extract the key alu elements, which with highest probability
  # for each position in the data frame
  # then we merge the data frame with the whole data.
  current_directory <- paste0(file_path, "/selection by chr cleaned")
  output_path <- paste0(file_path, "/highest probability")
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