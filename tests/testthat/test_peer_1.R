library(testthat)
library(ADVHEIDEL)
library(coda)

# Test 1: Single MCMC Chain with default test
test_that("advheidel.diag works with a single MCMC chain", {
  single_chain <- rnorm(1000, mean = 0, sd = 1)
  mcmc_single <- coda::mcmc(single_chain)
  result_single <- advheidel.diag(mcmc_single)

  expect_s3_class(result_single, "advHeidelDiag")
  expect_true(!is.null(result_single$results))  # Ensure there are results
  expect_equal(ncol(result_single$mcmc_chain), 1)  # One chain
})

# Test 2: Multiple MCMC Chains with default test
test_that("advheidel.diag works with multiple MCMC chains", {
  chain1 <- rnorm(1000, mean = 1, sd = 4)
  chain2 <- rnorm(1000, mean = 2, sd = 3)
  mcmc_multi <- cbind(chain1, chain2)
  mcmc_obj <- coda::mcmc(mcmc_multi)
  result_multi <- advheidel.diag(mcmc_obj)

  expect_s3_class(result_multi, "advHeidelDiag")
  expect_true(!is.null(result_multi$results))  # Ensure there are results
  expect_equal(ncol(result_multi$mcmc_chain), 2)  # Two chains
})

# Test 3: Single chain with CVM and Shapiro-Wilk tests
test_that("advheidel.diag works with CVM and Shapiro-Wilk tests", {
  single_chain <- rnorm(1000, mean = 0, sd = 1)
  mcmc_single <- coda::mcmc(single_chain)

  result_cvm <- advheidel.diag(mcmc_single, test_type = "cvm")
  result_shapiro <- advheidel.diag(mcmc_single, test_type = "shapiro")

  expect_s3_class(result_cvm, "advHeidelDiag")
  expect_s3_class(result_shapiro, "advHeidelDiag")
  expect_true(!is.null(result_cvm$results))  # Ensure there are results
  expect_true(!is.null(result_shapiro$results))  # Ensure there are results
})

# Test 4: Plot generation check (more robust)
test_that("advheidel.diag generates plots without errors", {
  chain1 <- rnorm(1000, mean = 1, sd = 4)
  chain2 <- rnorm(1000, mean = 2, sd = 3)
  mcmc_multi <- cbind(chain1, chain2)
  mcmc_obj <- coda::mcmc(mcmc_multi)
  result_multi <- advheidel.diag(mcmc_obj)

  # Ensure that the plot function does not throw errors
  expect_silent(plot(result_multi))
})
