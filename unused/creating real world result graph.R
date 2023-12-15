library(readxl)
library(dplyr)
library(stringr)
library(data.table)
library(igraph)
library(leiden)
library(RColorBrewer)
library(Rcpp)

# Load Hi-C data
setwd("your directory of data")

hic_data <- read.delim("hic_with_alu.txt", header = TRUE)

# Replace empty string with NA
# hic_data$unique_alu[hic_data$unique_alu == ""] <- NA
# Split the string at every comma and create the list
# list_alu <- strsplit(hic_data$unique_alu, split = ",")
# Assign list to a new column in hic_data
# hic_data$alu_list <- list_alu

hic_data$chromosome2 <- sub("^([^.]*)\\..*", "\\1", hic_data$frag2)

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

# you can check the number of non-zeros in the dataframe
# sum(hic_data$node_weights != 0)
vertices <- unique(c(hic_data$frag1, hic_data$frag2))

# we find the vertices are significantly low
edges <- data.frame(from = hic_data$frag1, 
                    to = hic_data$frag2, 
                    weight = hic_data$node_weights)
g <- graph_from_data_frame(d = edges, vertices = vertices, directed = FALSE)

V(g)$name <- vertices 
E(g)$weight <- edges$weight

resolution_param <- quantile(strength(g))[2] / (gorder(g) - 1) # This is a starting point, may need to adjust
cluster <- cluster_leiden(g, objective_function = "CPM", resolution_parameter = resolution_param, n_iterations = 100)

num_clusters <- length(unique(cluster$membership))
sizes <- table(cluster$membership)
size_of_largest_cluster <- max(sizes)

cat("Number of clusters:", num_clusters, "\n")
cat("Size of the largest cluster:", size_of_largest_cluster, "\n")

V(g)$community <- membership(cluster)

layout <- layout_with_fr(g)
V(g)$label <- NA

# Determine the predominant ALU family for each edge
hic_data$predominant_alu <- apply(hic_data[, c("AluJ_count", "AluS_count", "AluY_count")], 1, function(x) names(which.max(x)))

# Improved color scheme with high contrast and colorblind friendly
alu_colors <- c("AluJ_count" = "purple", "AluS_count" = "black", "AluY_count" = "red")

# Increase the scaling factor for edge width for better visibility
edge_scaling_factor <- 10

# Add alpha to edge colors for less visual clutter
edge_alpha <- 0.8
edges$color <- adjustcolor(alu_colors[hic_data$predominant_alu[match(paste(edges$from, edges$to), paste(hic_data$frag1, hic_data$frag2))]], alpha.f = edge_alpha)

# Use a layout that better separates the clusters, if available
# layout <- some_layout_function(g)  # Replace with a specific layout function if needed

# Plot the graph with updated parameters
png(filename = "high_res_plot.png", width = 2000, height = 2000, res = 300)

# Your existing plot code
plot(g,
     vertex.size = 4,
     vertex.label = NA,
     vertex.color = pal[V(g)$community],
     edge.width = E(g)$weight * 8,
     edge.color = edges$color,
     main = "Graph Cluster Plot",
     layout = layout
)

# Adding a legend (if your plot requires it)
legend("topright", legend = names(alu_colors), col = alu_colors, pch = 15)

# Turn off the PNG device, which saves the file
dev.off()

############ finish of the codes ############