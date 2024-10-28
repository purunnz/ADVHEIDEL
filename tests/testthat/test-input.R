library(testthat)
library(ADVHEIDEL)
library(coda)

set.seed(123)
mcmc_data <- cbind(rnorm(1000, mean = 1, sd = 4), rnorm(1000, mean = 2, sd = 3), rnorm(1000, mean = 3, sd = 2))
mcmc_obj <- coda::mcmc(mcmc_data)

# Test 1: Check if the function returns an object of class 'advHeidelDiag'
test_that("Function returns an object of class 'advHeidelDiag'", {
  result <- advheidel.diag(mcmc_obj, alpha = 0.05, eps = 0.1, test_type = "ks")
  expect_s3_class(result, "advHeidelDiag")
})

# Test 2: The 'results' field is a data frame
test_that("The 'results' field is a data frame", {
  result <- advheidel.diag(mcmc_obj, alpha = 0.05, eps = 0.1, test_type = "ks")
  expect_true(is.data.frame(result$results))
})

# Test 3: Number of rows in results matches the number of chains
test_that("Number of rows in results matches number of chains", {
  result <- advheidel.diag(mcmc_obj, alpha = 0.05, eps = 0.1, test_type = "ks")
  expect_equal(nrow(result$results), ncol(mcmc_obj))
})

# Test 4: Results contain expected columns
test_that("Results contain expected columns", {
  result <- advheidel.diag(mcmc_obj, alpha = 0.05, eps = 0.1, test_type = "ks")
  expected_columns <- c("Chain", "StationarityTest", "StartIteration", "PValue", "TestType", "HalfWidthTest", "Mean", "HalfWidth", "HalfWidthRatio")
  expect_true(all(expected_columns %in% colnames(result$results)))
})

# Test 5: Invalid input type
test_that("Input must be a matrix, data frame, or MCMC object", {
  expect_error(advheidel.diag(list(1, 2, 3)), "must be a matrix, data frame, or MCMC object")
})

# Test 6: Invalid stationarity test option triggers error
test_that("Invalid stationarity test option triggers error", {
  expect_error(advheidel.diag(mcmc_obj, test_type = "invalid_test"), "Unsupported test type. Choose 'ks', 'cvm', or 'shapiro'.")
})

# Test 7: Kolmogorov-Smirnov test runs correctly
test_that("Kolmogorov-Smirnov test runs correctly", {
  result <- advheidel.diag(mcmc_obj, test_type = "ks")
  expect_true("ks" %in% unique(result$results$TestType))
  expect_true(all(result$results$TestType == "ks"))
})

# Test 8: Cramer-von Mises test runs correctly
test_that("Cramer-von Mises test runs correctly", {
  result <- advheidel.diag(mcmc_obj, test_type = "cvm")
  expect_true("cvm" %in% unique(result$results$TestType))
  expect_true(all(result$results$TestType == "cvm"))
})

# Test 9: Shapiro-Wilk test runs correctly
test_that("Shapiro-Wilk test runs correctly", {
  result <- advheidel.diag(mcmc_obj, test_type = "shapiro")
  expect_true("shapiro" %in% unique(result$results$TestType))
  expect_true(all(result$results$TestType == "shapiro"))
})

# Test 10: Single chain case
test_that("Function works with a single chain", {
  single_chain <- rnorm(1000, mean = 0, sd = 1)
  mcmc_single <- coda::mcmc(single_chain)
  result_single <- advheidel.diag(mcmc_single, alpha = 0.05, eps = 0.1, test_type = "ks")
  expect_true(is.data.frame(result_single$results))
  expect_equal(nrow(result_single$results), 1)
})
