## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ADVHEIDEL)
library(goftest)

## -----------------------------------------------------------------------------
# Sample MCMC chains
set.seed(123)
chain1 <- rnorm(1000, mean = 1, sd = 4)
chain2 <- rnorm(1000, mean = 2, sd = 3)
mcmc_data <- cbind(chain1, chain2)

## -----------------------------------------------------------------------------

# Apply the Heidelberger-Welch diagnostic
result <- advheidel.diag(mcmc_data, custommetrics = TRUE)
print(result)
summary(result)

## -----------------------------------------------------------------------------
# Plot diagnostic results
plot(result)

## ----eval=FALSE---------------------------------------------------------------
#  # Example test case
#  test_that("heidelWelchMetric calculates mean, variance, and half-width correctly", {
#    x <- rnorm(100, mean = 5, sd = 2)
#    result <- heidelWelchMetric(x)
#    expect_equal(result$Mean[1], mean(x), tolerance = 1e-8)
#  })

