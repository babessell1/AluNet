library(testthat)
library(Rcpp)

sourceCpp('compile_test.cpp')  

test_that("compile_success", {
  nodes <- list("Node1", "Node2", "Node3")
  edges <- list(c("Node1", "Node2"), c("Node2", "Node3"))
  weights <- c(1.5, 2.5)
  
  graph <- compile_success(nodes, edges, weights)
  
  expect_true(is.list(graph))
  expect_equal(graph$nodes, nodes)
  expect_equal(graph$edges, edges)
  expect_equal(graph$weights, weights)
  
  expect_equal(length(graph$nodes), length(nodes))
  expect_equal(length(graph$edges), length(edges))
  expect_equal(length(graph$weights), length(weights))
})
