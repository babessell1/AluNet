library(dplyr)
library(purrr)
library(readr)
library(data.table)
library(stringr)
library(readxl)
library(RColorBrewer)

# change your working directory
# setwd()
current_directory <- getwd()

# place all your data to the uncleaned data folder
file_path <- paste0(current_directory, "/uncleaned data")
output_path <- paste0(current_directory, "/selection by chr")
# Create the directory
if (!file.exists(output_path)) {
  dir.create(output_path)
}

setwd(file_path)

hic_data <- read.delim("hic.txt", header = TRUE)
hic_data$chromosome <- sub("^([^.]*)\\..*", "\\1", hic_data$frag1)
unique_chromosomes <- unique(hic_data$chromosome)

read_and_filter <- function(file, chromosome) {
  data <- read.csv(file)
  filtered_data <- data %>% filter(X.seq_name == chromosome)
  return(filtered_data)
}

files <- list.files(file_path, pattern="*.csv", full.names = TRUE)

# we first extract the chromosome's alu elements

for(i in unique_chromosomes){
  cat(paste0("now we are working on chromosome: ", i))
  
  merged_data <- map_df(files,  ~read_and_filter(., i))
  
  print(nrow(merged_data))
  
  output_file <- paste0(output_path, "/", i)
  write.csv(merged_data, output_file)
}

# we first extract the key alu elements, which with highest probability
# for each position in the data frame
# then we merge the data frame with the whole data.
file_path <- output_path
output_path <- paste0(current_directory, "/highest probability")
# Create the directory
if (!file.exists(output_path)) {
  dir.create(output_path)
}

setwd(file_path)

files <- list.files(file_path, pattern="*", full.names = TRUE)

for(i in files){
  last_item <- basename(i)
  last_item
  
  cat(paste0("now we are working on chromosome: ", last_item))
  
  test_file <- read.csv(i)
  test_file <- data.table(test_file)
  
  df_selected <- test_file[, .SD[which.min(e.value)], by = .(ali.st, ali.en)]
  print(nrow(df_selected))
  write.csv(df_selected, paste0(output_path, '/', last_item))
}

# finally, we merge the alu elements with the hic data
# the hic data frame contains the 

file_path <- paste0(current_directory, "/uncleaned data")

setwd(file_path)

# Load Hi-C data
hic_data <- read.delim("hic.txt", header = TRUE)
matches <- str_match(hic_data$frag1, "chr[^.]+\\.([0-9]+)")
matches2 <- str_match(hic_data$frag2, "chr[^.]+\\.([0-9]+)")
# Extract the captured numeric part and convert to numeric
numbers <- as.numeric(matches[, 2])
numbers2 <- as.numeric(matches2[, 2])
# Apply the floor operation in a vectorized manner
hic_data$start <- ifelse(!is.na(numbers), floor(numbers / 5000) * 5000, NA)
hic_data$start2 <- ifelse(!is.na(numbers2), floor(numbers2 / 5000) * 5000, NA)
# Check the number of NA values
sum(is.na(hic_data$start))
hic_data$end <- hic_data$start + 5000
sum(is.na(hic_data$start2))
hic_data$end2 <- hic_data$start2 + 5000

hic_data$chromosome <- sub("^([^.]*)\\..*", "\\1", hic_data$frag1)
unique_chromosomes <- unique(hic_data$chromosome)
print(unique_chromosomes)

file_path <- paste0(current_directory, "/highest probability")

files <- list.files(file_path, pattern="*", full.names = TRUE)

setwd(file_path)

big_df <- data.table()
for(file in files){
  chromo <- basename(file)
  cat(paste0("now we are working on chromosome: ", chromo), "\n")
  
  dt <- fread(file)
  dt[, 'chromo' := chromo] # Create a new column 'chromo' with the filename
  
  # Bind the new data table to the 'big_df'
  big_df <- rbind(big_df, dt)
}

big_df <- rename(big_df, end = ali.en, start = ali.st)

# set table ready to merge
setDT(hic_data)
setDT(big_df)

## now convert the family of ALU to subfamilies

unique_alu_types <- unique(big_df$family_name)
print(length(unique_alu_types))

big_df$family_name <- gsub("^AluJ.*", "AluJ", big_df$family_name) 
big_df$family_name <- gsub("^AluS.*", "AluS", big_df$family_name) 
big_df$family_name <- gsub("^AluY.*", "AluY", big_df$family_name)

unique_alu_types_revised <- unique(big_df$family_name)

unique_alu_types_revised

# generate a index of hic data and using the hic data to operate the merging

hic_data[, index := .I]

overlap_frag1 <- big_df[hic_data, on = .(chromo = chromosome, start >= start, end <= end),
                        .(mapping = i.index, frag1 = i.frag1, alu_name = x.family_name),
                        nomatch = 0L, allow.cartesian = TRUE]

overlap_frag2 <- big_df[hic_data, on = .(chromo = chromosome, start >= start2, end <= end2),
                        .(mapping = i.index, frag2 = i.frag2, alu_name = x.family_name),
                        nomatch = 0L, allow.cartesian = TRUE]

hic_data[, c("AluJ_count", "AluS_count", "AluY_count", "alu_total") := 0]

overlap_frag1_agg <- overlap_frag1[, .(aluJ_count = sum(alu_name=="AluJ"), 
                                       aluS_count = sum(alu_name=="AluS"), 
                                       aluY_count = sum(alu_name=="AluY")), 
                                   by = mapping]

overlap_frag2_agg <- overlap_frag2[, .(aluJ_count = sum(alu_name=="AluJ"), 
                                       aluS_count = sum(alu_name=="AluS"), 
                                       aluY_count = sum(alu_name=="AluY")), 
                                   by = mapping]

merged_overlap <- merge(overlap_frag1_agg, overlap_frag2_agg, by = "mapping", all = TRUE)
merged_overlap[is.na(merged_overlap)] <- 0 # replace NAs in the merged dataset with 0s
merged_overlap[, AluJ_count := pmin(aluJ_count.x, aluJ_count.y)]
merged_overlap[, AluS_count := pmin(aluS_count.x, aluS_count.y)]
merged_overlap[, AluY_count := pmin(aluY_count.x, aluY_count.y)]
merged_overlap[, alu_total := AluJ_count + AluS_count + AluY_count]


count_alu_subset <- merged_overlap[, .(AluJ_count = AluJ_count, 
                                       AluS_count = AluS_count,
                                       AluY_count = AluY_count,
                                       alu_total = alu_total), by = mapping]

setnames(count_alu_subset, old = "mapping", new = "index")
hic_data[count_alu_subset, on = "index", `:=` (AluJ_count = i.AluJ_count, 
                                               AluS_count = i.AluS_count, 
                                               AluY_count = i.AluY_count, 
                                               alu_total = i.alu_total)]

intermediate_result <- paste0(current_directory, "/cleaned hic with alu")
# Create the directory
if (!file.exists(intermediate_result)) {
  dir.create(intermediate_result)
}

setwd(intermediate_result)
write.table(hic_data, file = "hic_with_alu.txt", sep = "\t")

setwd(current_directory)
rm(list = ls())

file_path <- paste0(current_directory, "/cleaned hic with alu")
setwd(file_path)

hic_data <- read.delim("hic_with_alu.txt", header = TRUE) 
hic_data$chromosome2 <- sub("^([^.]*)\\..*", "\\1", hic_data$frag2)

file_path <- paste0(current_directory, "/uncleaned data")
setwd(file_path)

# Load ATAC-seq data
atac_data <- read_excel("atac.xlsx",col_names = TRUE)
atac_data <- atac_data[-(1:18),]
colnames(atac_data) <- c("Peak_ID", "hg38_Chromosome", "start", "end", 
                         "Peak_Width", "Annotation", "Distance_To_TSS", 
                         "Gene_Symbol", "GC_Percent", "CTCF")
unique_chromosomes <- unique(atac_data$hg38_Chromosome)
atac_data$start <- as.numeric(atac_data$start)
atac_data$end <- as.numeric(atac_data$end)

setDT(hic_data)
setDT(atac_data)

# for chromosome, frag1
hic_data[, atac_peak := 0]  # Initialize the column with zeros
overlap <- atac_data[hic_data, on = .(hg38_Chromosome = chromosome, start >= start, end <= end), nomatch = 0L, allow.cartesian = TRUE,
                     .(frag1 = i.frag1, start = i.start, end = i.end, atac_peak = 1)]
overlap[, atac_peak := 1]
hic_data[overlap, atac_peak := i.atac_peak, on = .(frag1, start, end)]
nrow(hic_data[atac_peak == 1])

# for chromosome, frag2
hic_data[, atac_peak2 := 0]  # Initialize the column with zeros
overlap2 <- atac_data[hic_data, on = .(hg38_Chromosome = chromosome2, start >= start2, end <= end2), nomatch = 0L, allow.cartesian = TRUE,
                      .(frag2 = i.frag2, start2 = i.start2, end2 = i.end2, atac_peak2 = 1)]
overlap2[, atac_peak2 := 1]
hic_data[overlap2, atac_peak2 := i.atac_peak2, on = .(frag2, start2, end2)]
nrow(hic_data[atac_peak2 == 1])

# combine together
hic_data$peak <- (hic_data$atac_peak+hic_data$atac_peak2)
hic_data$peak <- ifelse(hic_data$peak == 2, 1, 0)
hic_data[, c("atac_peak", "atac_peak2") := NULL]
setnames(hic_data, "peak", "atac_peak")

# combine together
setnames(hic_data, "alu_total", "alu_peak")

# ##########################
# revisions
hic_data$node_weights <- hic_data$atac_peak * hic_data$rawCount * hic_data$alu_peak

hic_data <- hic_data[which(hic_data$chromosome == "chr22")]

# print(summary(hic_data$node_weights))
## we first try to give each weight a very small value to let every data in the model
# hic_data$node_weights <- hic_data$node_weights + 0.01
hic_data <- hic_data[which(hic_data$node_weights != 0), ]
hic_data$node_weights <- hic_data$node_weights / max(hic_data$node_weights)

final_result <- paste0(current_directory, "/final output")
# Create the directory
if (!file.exists(final_result)) {
  dir.create(final_result)
}
setwd(final_result)
write.csv(hic_data, file = "edges_data_frame.csv", row.names = TRUE)