library(testthat)
library(Rcpp)

sourceCpp('add_test.cpp') 

test_that("addTwoIntegers", {
  expect_equal(addTwoIntegers(1, 2), 3)
  expect_equal(addTwoIntegers(-1, 1), 0)
  expect_equal(addTwoIntegers(-1, -2), -3)
})

