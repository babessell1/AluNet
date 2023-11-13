library(Rcpp)
library(igraph)
library(testthat)
library(leidenAlg)

sourceCpp("../rcpp/Leiden.cpp") # ues the rcpp file that is located in another directory

generate_large_weighted_graph <- function(n, p) {
  graph <- erdos.renyi.game(n, p, type = "gnp", directed = FALSE)
  E(graph)$weight <- runif(ecount(graph))  
  return(graph)
}

igraph_to_list <- function(graph) {
  edges <- get.edgelist(graph)
  weights <- E(graph)$weight
  if (!is.null(weights)) {
    edge_list <- cbind(edges, weights)
  } else {
    edge_list <- edges
  }
  graph_list <- split(edge_list, rep(1:nrow(edge_list), each = ncol(edge_list)))
  return(graph_list)
}

create_community <- function(index, nodes) {
  community <- Rcpp::new(Community, as.integer(nodes), as.integer(index))
  return(community)
}

create_partition <- function(communities) {
  partition <- Rcpp::new(Partition, communities)
  return(partition)
}

# Generate a large weighted graph using an adjacency matrix
generate_large_weighted_graph_adj_matrix <- function(n, p) {
  adj_matrix <- matrix(runif(n * n) < p, n, n)
  diag(adj_matrix) <- 0
  adj_matrix[lower.tri(adj_matrix)] <- t(adj_matrix)[lower.tri(adj_matrix)]
  adj_matrix[adj_matrix == TRUE] <- runif(sum(adj_matrix == TRUE))
  graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)
  return(graph)
}

# Convert igraph object to adjacency list
igraph_to_list_adj <- function(graph) {
  adj_list <- as_adj_list(graph, mode = "all", attr = "weight")
  return(adj_list)
}

# Create community partition using Leiden algorithm
create_community_partition_leiden <- function(graph) {
  community <- leiden(graph)
  return(community)
}

test_that("Community::aggregateWeights works correctly", {
  large_graph <- generate_large_weighted_graph(10000, 0.01)
  community <- create_community(1, sample(1:10000, 500)) 
  weight_sum <- community$aggregateWeights(large_graph)
  expect_true(is.numeric(weight_sum))
  
  total_weight <- sum(E(large_graph)$weight)
  expect_true(weight_sum <= total_weight)
})

test_that("Community::size returns correct size", {
  community <- create_community(1, 1:500) # Community with 500 nodes
  expect_equal(community$size(), 500)
})

test_that("Partition constructor works correctly", {
  communities <- list(create_community(1, 1:500), create_community(2, 501:1000))
  partition <- create_partition(communities)
  expect_true(length(partition$getCommunityIndices()) == 2)
  
  expect_true(all(sapply(partition$communityIndexMap, function(x) x$communityIndex %in% c(1, 2))))
})

test_that("Partition::getCommunityIndices works correctly", {
  communities <- list(create_community(1, 1:500), create_community(2, 501:1000))
  partition <- create_partition(communities)
  expect_equal(partition$getCommunityIndices(), c(1, 2))
})

test_that("Partition::flattenPartition works correctly", {
  communities <- list(create_community(1, 1:500), create_community(2, 501:1000))
  partition <- create_partition(communities)
  partition$flattenPartition()
  
  for (community_index in partition$getCommunityIndices()) {
    expect_false(any(class(partition$communityIndexMap[[community_index]]) == "Community"))
  }
})

test_that("Partition::addCommunity adds a community correctly", {
  partition <- create_partition(list())
  community <- create_community(3, 1001:1500)
  partition$addCommunity(community)
  expect_true(3 %in% partition$getCommunityIndices())
  
  expect_true(all(1001:1500 %in% partition$communityIndexMap[[3]]$nodeIndices))
})

test_that("Partition::updateCommunityMembershipSearch updates correctly", {
  communities <- list(create_community(1, 1:500), create_community(2, 501:1000))
  partition <- create_partition(communities)
  test_node <- 250
  partition$updateCommunityMembershipSearch(test_node, 2)
  
  expect_false(test_node %in% partition$communityIndexMap[[1]]$nodeIndices)
  expect_true(test_node %in% partition$communityIndexMap[[2]]$nodeIndices)
})

test_that("Partition::updateCommunityMembership updates correctly", {
  communities <- list(create_community(1, 1:500), create_community(2, 501:1000))
  partition <- create_partition(communities)
  test_node <- 250
  partition$updateCommunityMembership(250, 1, 2)
  # Check for correct membership update
  expect_false(test_node %in% partition$communityIndexMap[[1]]$nodeIndices)
  expect_true(test_node %in% partition$communityIndexMap[[2]]$nodeIndices)
})

test_that("Leiden algorithm results are comparable to iGraph's", {
  small_graph <- make_graph("Zachary") 
  communities_leiden <- runLeiden(igraph_to_list(small_graph), iterations = 10)
  igraph_communities <- cluster_louvain(small_graph)
  # Compare community structures
  expect_true(all(sort(communities_leiden$communities) == sort(igraph_communities$membership)))
})

test_that("Optimizer class functions work correctly", {
  large_graph <- generate_large_weighted_graph(10000, 0.01)
  communities <- list(create_community(1, 1:500), create_community(2, 501:1000))
  partition <- create_partition(communities)
  gamma <- 1.0  
  optimizer <- Rcpp::new(Optimizer, large_graph, partition, gamma)
  
  expect_true(optimizer$moveNodesFast())
  
  refined_partition <- optimizer$refinePartition()
  expect_true(is(refined_partition, "Partition"))  
  
  subset <- create_community(3, 1001:1500)  
  merged_partition <- optimizer$mergeNodesSubset(subset)
  expect_true(is(merged_partition, "Partition"))  # Check if merged_partition is a valid Partition object
  
  aggregated_graph <- optimizer$aggregateGraph(large_graph, partition)
  expect_true(is(aggregated_graph, "Graph"))  # Check if aggregated_graph is a valid Graph object
  
  potts_value <- optimizer$constantPotts(gamma)
  expect_true(is.numeric(potts_value))
  
  optimizer$optimize()
  final_partition <- optimizer$P
  expect_true(length(final_partition$getCommunityIndices()) > 1)
})

test_that("Generated graph from adjacency matrix has correct properties", {
  graph <- generate_large_weighted_graph_adj_matrix(100, 0.1)
  expect_true(is_undirected(graph))
  expect_true(is_weighted(graph))
  expect_true(gorder(graph) == 100)
})

test_that("Leiden algorithm detects communities from adjacency matrix graph", {
  graph <- generate_large_weighted_graph_adj_matrix(100, 0.1)
  partition <- create_community_partition_leiden(graph)
  expect_true(length(partition$membership) == gorder(graph))
})

test_that("Adjacency list conversion maintains structure for adjacency matrix graph", {
  graph <- generate_large_weighted_graph_adj_matrix(100, 0.1)
  adj_list <- igraph_to_list_adj(graph)
  expect_true(length(adj_list) == gorder(graph))
})

test_that("Leiden algorithm results from adjacency matrix graph are comparable to iGraph's Louvain", {
  graph <- generate_large_weighted_graph_adj_matrix(100, 0.1)
  communities_leiden <- create_community_partition_leiden(graph)
  communities_louvain <- cluster_louvain(graph)
  expect_true(all(sort(communities_leiden$membership) == sort(communities_louvain$membership)))
})



