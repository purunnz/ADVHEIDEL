library(testthat)
library(ADVHEIDEL)
library(coda)

# Test 1: Testing with a data frame input
test_that("advheidel.diag works with data.frame input", {
  mcmc_data_df <- data.frame(chain1 = rnorm(1000, mean = 1, sd = 2), chain2 = rnorm(1000, mean = 2, sd = 1))
  result_df_input <- advheidel.diag(mcmc_data_df)
  expect_s3_class(result_df_input, "advHeidelDiag")
  expect_true(!is.null(result_df_input$results))
})

# Test 2: Testing with matrix input
test_that("advheidel.diag works with matrix input", {
  mcmc_data_matrix <- matrix(rnorm(2000), ncol = 2)
  result_matrix_input <- advheidel.diag(mcmc_data_matrix)
  expect_s3_class(result_matrix_input, "advHeidelDiag")
  expect_true(!is.null(result_matrix_input$results))
})

# Test 3: Testing with a single vector
test_that("advheidel.diag works with single vector", {
  mcmc_chain <- rnorm(1000, mean = 0, sd = 1)
  result_single_vector <- advheidel.diag(mcmc_chain)
  expect_s3_class(result_single_vector, "advHeidelDiag")
  expect_true(!is.null(result_single_vector$results))
})

# Test 4: Handling invalid test type
test_that("advheidel.diag handles invalid test type", {
  mcmc_single <- coda::mcmc(rnorm(1000))
  expect_error(advheidel.diag(mcmc_single, test_type = "invalid_test"), "Unsupported test type. Choose 'ks', 'cvm', or 'shapiro'.")
})

# Test 5: Handling missing chain input
test_that("advheidel.diag handles missing chain input gracefully", {
  expect_error(advheidel.diag(), "argument \"mcmc_chain\" is missing")
})

# Test 6: Handling NA values in the MCMC chains
test_that("advheidel.diag stops with an error when NA values are present in the chain", {
  mcmc_chain_with_na <- matrix(rnorm(1000), ncol = 1)
  mcmc_chain_with_na[100] <- NA  # Insert an NA value
  expect_error(advheidel.diag(mcmc_chain_with_na), "The 'mcmc_chain' contains Inf, NA, or NaN values. Please provide valid data.")
})

# Test 7: Handling Inf values in the MCMC chains
test_that("advheidel.diag stops with an error when Inf values are present in the chain", {
  chain_with_inf <- rnorm(1000)
  chain_with_inf[c(100, 200, 500)] <- Inf
  chain_with_na2 <- rnorm(1000)
  mcmc_with_inf <- cbind(chain_with_inf, chain_with_na2)
  mcmc_with_inf <- coda::mcmc(mcmc_with_inf)
  expect_error(advheidel.diag(mcmc_with_inf), "The 'mcmc_chain' contains Inf, NA, or NaN values. Please provide valid data.")
})

# Test 8: Handling Inf, NA, or NaN values in the chain
test_that("advheidel.diag handles Inf, NA, or NaN values in the chain", {
  mcmc_with_na <- matrix(c(1, 2, NA, 4, Inf, NaN), ncol = 1)
  expect_error(advheidel.diag(mcmc_with_na), "The 'mcmc_chain' contains Inf, NA, or NaN values. Please provide valid data.")
})

# Test 9: Testing dependent chains
test_that("advheidel.diag works with dependent chains", {
  chain_a <- rnorm(1000)
  chain_b <- 2 * chain_a + rnorm(1000, 0, 0.1)  # Dependent on chain_a
  mcmc_dependent <- cbind(chain_a, chain_b)
  result_dependent <- advheidel.diag(mcmc_dependent)
  expect_s3_class(result_dependent, "advHeidelDiag")
  expect_true(!is.null(result_dependent$results))
})

# Test 10: Performance test with a large MCMC chain
test_that("advheidel.diag works with a large MCMC chain", {
  large_chain <- rnorm(100000)  # 100,000 data points
  mcmc_large <- coda::mcmc(large_chain)

  time_taken <- system.time({
    result_large <- advheidel.diag(mcmc_large)
  })

  print(time_taken)  # Print the time taken for execution

  expect_s3_class(result_large, "advHeidelDiag")
  expect_true(!is.null(result_large$results))
})
