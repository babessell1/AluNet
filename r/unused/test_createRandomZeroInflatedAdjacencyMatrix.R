
library(testthat)

source('createRandomZeroInflatedAdjacencyMatrix.r')

test_that("Function returns a matrix of correct dimensions", {
  num_vertices <- 5
  result <- createRandomZeroInflatedAdjacencyMatrix(num_vertices, 0.3)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(num_vertices, num_vertices))
})

test_that("Matrix has correct proportion of zeros", {
  num_vertices <- 5
  zero_prob <- 0.3
  result <- createRandomZeroInflatedAdjacencyMatrix(num_vertices, zero_prob)
  zero_count <- sum(result == 0)
  total_elements <- num_vertices * num_vertices
  observed_zero_prob <- zero_count / total_elements
  expect_equal(observed_zero_prob, zero_prob, tolerance = 0.1)
})

test_that("Adjacency matrix is symmetric", {
  num_vertices <- 5
  result <- createRandomZeroInflatedAdjacencyMatrix(num_vertices, 0.3)
  expect_true(identical(result, t(result)))
})

test_that("Non-zero weights are positive", {
  num_vertices <- 5
  result <- createRandomZeroInflatedAdjacencyMatrix(num_vertices, 0.3)
  non_zero_weights <- result[result != 0]
  expect_true(all(non_zero_weights > 0))
})
