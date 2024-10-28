library(testthat)
library(ADVHEIDEL)
library(Rcpp)
library(goftest)

# sourceCpp("src/heidelWelchMetric.cpp")


# Test 1: Basic functionality with a small numeric vector
test_that("heidelWelchMetric works with a basic vector", {
  # Sample data for testing
  x <- c(5, 6, 7, 8, 9)
  result <- heidelWelchMetric(x)

  # Loop through each discard level result in heidelWelchMetric output
  for (i in seq_along(result$Discard.Ratio)) {
    # Extract expected values
    discard_ratio <- result$Discard.Ratio[i]
    trimmed_x <- x[(floor(length(x) * discard_ratio) + 1):length(x)]

    # Calculate expected values based on trimmed data
    expected_mean <- mean(trimmed_x)
    expected_variance <- var(trimmed_x)
    expected_halfwidth <- qt(1 - 0.05 / 2, df = length(trimmed_x) - 1) * sqrt(expected_variance / length(trimmed_x))

    # Check mean and variance
    expect_equal(result$Mean[i], expected_mean, tolerance = 1e-8)
    expect_equal(result$Variance[i], expected_variance, tolerance = 1e-8)

    # Check half-width and p-value with adjusted tolerance for minor differences
    expect_equal(result$Half.Width[i], expected_halfwidth, tolerance = 1e-8)
    expect_equal(result$p.value[i], 1, tolerance = 0.01)  # Adjusted tolerance to 0.01 for p-value

    # Check half-width ratio and stationarity results based on predefined thresholds
    expect_equal(result$Half.Width.Ratio[i], abs(expected_halfwidth / expected_mean), tolerance = 1e-8)
    expect_true(result$Stationarity.Test[i])  # Assuming test passes for this example data
  }
})

# Test 2: Functionality with a normally distributed vector
test_that("heidelWelchMetric calculates mean, variance, and half-width correctly", {
  set.seed(123)
  x <- rnorm(100, mean = 5, sd = 2)
  result <- heidelWelchMetric(x)

  # Loop through each discard level to validate expected mean, variance, and half-width
  for (i in seq_along(result$Discard.Ratio)) {
    discard_ratio <- result$Discard.Ratio[i]
    trimmed_x <- x[(floor(length(x) * discard_ratio) + 1):length(x)]

    expected_mean <- mean(trimmed_x)
    expected_variance <- var(trimmed_x)
    expected_halfwidth <- qt(1 - 0.05 / 2, df = length(trimmed_x) - 1) * sqrt(expected_variance / length(trimmed_x))

    expect_equal(result$Mean[i], expected_mean, tolerance = 1e-8)
    expect_equal(result$Variance[i], expected_variance, tolerance = 1e-8)
    expect_equal(result$Half.Width[i], expected_halfwidth, tolerance = 1e-8)
  }
})

# Test 3: Single value vector handling
test_that("heidelWelchMetric handles single value vector", {
  x <- c(10)
  result <- heidelWelchMetric(x)

  # Ensure correct handling of single value with expected NA values
  expect_equal(result$Mean[1], 10)
  expect_true(is.na(result$Variance[1]))
  expect_true(is.na(result$Half.Width[1]))
})

# Test 4: Empty vector handling
test_that("heidelWelchMetric raises an error for empty vector", {
  x <- numeric(0)
  expect_error(heidelWelchMetric(x), "Input vector must contain at least one value.")
})

# Test 5: Custom alpha level and its effect on half-width
test_that("heidelWelchMetric adjusts half-width for different alpha levels", {
  set.seed(123)
  x <- rnorm(100, mean = 5, sd = 2)

  result_default <- heidelWelchMetric(x, alpha = 0.05)
  result_custom <- heidelWelchMetric(x, alpha = 0.01)

  # Validate half-width changes for custom alpha level across discard ratios
  for (i in seq_along(result_default$Discard.Ratio)) {
    discard_ratio <- result_default$Discard.Ratio[i]
    trimmed_x <- x[(floor(length(x) * discard_ratio) + 1):length(x)]

    # Expected half-width calculations for each alpha level
    expected_halfwidth_default <- qt(1 - 0.05 / 2, df = length(trimmed_x) - 1) * sqrt(var(trimmed_x) / length(trimmed_x))
    expected_halfwidth_custom <- qt(1 - 0.01 / 2, df = length(trimmed_x) - 1) * sqrt(var(trimmed_x) / length(trimmed_x))

    expect_equal(result_default$Half.Width[i], expected_halfwidth_default, tolerance = 1e-8)
    expect_equal(result_custom$Half.Width[i], expected_halfwidth_custom, tolerance = 1e-8)
  }
})

# Test 6: Handling large vectors to ensure stability
test_that("heidelWelchMetric handles large vectors", {
  set.seed(123)
  x <- rnorm(1e6, mean = 5, sd = 2)
  result <- heidelWelchMetric(x)

  # Validate against mean, variance, and half-width for large vector
  for (i in seq_along(result$Discard.Ratio)) {
    discard_ratio <- result$Discard.Ratio[i]
    trimmed_x <- x[(floor(length(x) * discard_ratio) + 1):length(x)]

    expected_mean <- mean(trimmed_x)
    expected_variance <- var(trimmed_x)
    expected_halfwidth <- qt(1 - 0.05 / 2, df = length(trimmed_x) - 1) * sqrt(expected_variance / length(trimmed_x))

    expect_equal(result$Mean[i], expected_mean, tolerance = 1e-8)
    expect_equal(result$Variance[i], expected_variance, tolerance = 1e-8)
    expect_equal(result$Half.Width[i], expected_halfwidth, tolerance = 1e-8)
  }
})

# Test 7: Custom metrics functionality in advheidel.diag
test_that("advheidel.diag works with custom metrics", {
  chain1 <- rnorm(1000)
  chain2 <- rnorm(1000)
  custom_chain <- cbind(chain1, chain2)
  mcmc_custom <- coda::mcmc(custom_chain)
  result_custom <- advheidel.diag(mcmc_custom, custommetrics = TRUE)

  # Validate class and presence of custom metrics data
  expect_s3_class(result_custom, "advHeidelDiag")
  expect_true(!is.null(result_custom$results))

  # Confirm presence of custom metrics by checking one chain's custom results
  expect_true("CustomMetric" %in% names(result_custom$results))
  expect_true(!is.null(result_custom$results$CustomMetric[[1]]))
})

