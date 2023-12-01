library(readxl)
library(dplyr)
library(stringr)
library(data.table)

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

hic_data$chromosome2 <- sub("^([^.]*)\\..*", "\\1", hic_data$frag2)
unique_chromosomes2 <- unique(hic_data$chromosome2)
print(unique_chromosomes2)





# Load ATAC-seq data
atac_data <- read_excel("atac.xlsx",col_names = TRUE)
atac_data <- atac_data[-(1:18),]
colnames(atac_data) <- c("Peak_ID", "hg38_Chromosome", "start", "end", 
                         "Peak_Width", "Annotation", "Distance_To_TSS", 
                         "Gene_Symbol", "GC_Percent", "CTCF")
unique_chromosomes <- unique(atac_data$hg38_Chromosome)
print(unique_chromosomes)
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

hic_data[, hic_order := .I]
atac_data[, atac_order := .I]

overlap <- atac_data[hic_data, 
                     on = .(hg38_Chromosome = chromosome, start >= start, end <= end), 
                     nomatch = 0L, 
                     allow.cartesian = TRUE,
                     .(hic_start = i.start, hic_end = i.end, hic_order = i.hic_order, 
                       atac_start = x.start, atac_end = x.end, atac_order = x.atac_order)]
unique(overlap$hic_order)

hic_data[, hic_order2 := .I]
overlap2 <- atac_data[hic_data, 
                     on = .(hg38_Chromosome = chromosome2, start >= start2, end <= end2), 
                     nomatch = 0L, 
                     allow.cartesian = TRUE,
                     .(hic_start2 = i.start2, hic_end2 = i.end2, hic_order = i.hic_order, 
                       atac_start = x.start, atac_end = x.end, atac_order = x.atac_order)]
unique(overlap2$hic_order)

# combine together
hic_data$peak <- pmax(hic_data$atac_peak, hic_data$atac_peak2)
hic_data[, c("atac_peak", "atac_peak2") := NULL]
setnames(hic_data, "peak", "atac_peak")






# Load ALU element data
alu_data <- read_excel("alu.xlsx",col_names = TRUE)
alu_data <- alu_data[-c(1:2),]
colnames(alu_data) <- c("Chromosome", "start", "end", "Strand", 
                         "Annotation", "Subfamily", "Nearest gene")

setDT(alu_data)
alu_data$start <- as.numeric(alu_data$start)
alu_data$end <- as.numeric(alu_data$end)

# for chromosome, frag1
hic_data[, alu_peak := 0]  # Initialize the column with zeros
overlap_alu <- alu_data[hic_data, on = .(Chromosome = chromosome, start >= start, end <= end), nomatch = 0L, allow.cartesian = TRUE]
overlap_alu[, alu_peak := 1]
hic_data[overlap_alu, alu_peak := i.alu_peak, on = .(frag1, start, end)]
nrow(hic_data[alu_peak == 1])

# for chromosome, frag2
hic_data[, alu_peak2 := 0]  # Initialize the column with zeros
overlap_alu2 <- alu_data[hic_data, on = .(Chromosome = chromosome2, start >= start2, end <= end2), nomatch = 0L, allow.cartesian = TRUE,
                         .(frag2 = i.frag2, start2 = i.start2, end2 = i.end2, atac_peak2 = 1)]
overlap_alu2[, alu_peak2 := 1]
hic_data[overlap_alu2, alu_peak2 := i.alu_peak2, on = .(frag2, start2, end2)]
nrow(hic_data[alu_peak2 == 1])


alu_data[, alu_order := .I]

overlap_alu <- alu_data[hic_data, 
                     on = .(Chromosome = chromosome, start >= start, end <= end), 
                     nomatch = 0L, 
                     allow.cartesian = TRUE,
                     .(hic_start = i.start, hic_end = i.end, hic_order = i.hic_order, 
                       alu_start = x.start, alu_end = x.end, alu_order = x.alu_order)]
unique(overlap_alu$hic_order)


overlap_alu2 <- alu_data[hic_data, 
                        on = .(Chromosome = chromosome2, start >= start2, end <= end2), 
                        nomatch = 0L, 
                        allow.cartesian = TRUE,
                        .(hic_start2 = i.start2, hic_end2 = i.end2, hic_order = i.hic_order, 
                          alu_start = x.start, alu_end = x.end, alu_order = x.alu_order)]
unique(overlap_alu2$hic_order)


# combine together
hic_data$peak <- pmax(hic_data$alu_peak, hic_data$alu_peak2)
hic_data[, c("hic_order", "hic_order2", "alu_peak", "alu_peak2") := NULL]
setnames(hic_data, "peak", "alu_peak")




# combine together

hic_data[, weight := atac_peak * alu_peak]