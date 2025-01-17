% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/advheidel.R
\name{advheidel.diag}
\alias{advheidel.diag}
\title{advheidel.diag Function: Enhanced Heidelberger-Welch MCMC Diagnostics}
\usage{
advheidel.diag(
  mcmc_chain,
  alpha = 0.05,
  eps = 0.1,
  test_type = "ks",
  custommetrics = FALSE
)
}
\arguments{
\item{mcmc_chain}{An MCMC object (matrix or data frame) containing one or more chains.}

\item{alpha}{Numeric. The significance level for the stationarity and half-width tests. Default is 0.05.}

\item{eps}{Numeric. The threshold for the Half-Width ratio test. Default is 0.1.}

\item{test_type}{Character. The stationarity test to apply. Options are "ks" for Kolmogorov-Smirnov,
"cvm" for Cramer-von Mises, "shapiro" for Shapiro-Wilk test. Default is "ks".}

\item{custommetrics}{Logical. If TRUE, the function calculates additional custom metrics using the
\code{heidelWelchMetric} function. Default is FALSE.}
}
\value{
A list object of class \code{advHeidelDiag} containing:
\describe{
  \item{results}{A data frame with diagnostic results for each chain.}
  \item{mcmc_chain}{The input MCMC chain.}
}
}
\description{
This function performs an advanced Heidelberger-Welch MCMC diagnostic on one or more chains
with options for different stationarity tests. It includes the Kolmogorov-Smirnov (KS), Cramer-von Mises (CVM),
Shapiro-Wilk tests, and custom metrics when applicable. The function also calculates the Half-Width ratio for each chain.
}
\details{
This function applies the Heidelberger-Welch stationarity and Half-Width tests to diagnose convergence issues in MCMC chains.
It allows users to choose from three different stationarity tests (\code{ks}, \code{cvm}, \code{shapiro}) and provides additional metrics such as the Half-Width ratio.
The \code{custommetrics} parameter, if set to TRUE, enables the calculation of additional diagnostic metrics using the C++ implementation via the \code{heidelWelchMetric} function.

The stationarity test involves discarding parts of the MCMC chains and checking for convergence. If custom metrics are enabled, the function returns diagnostic results, including mean, variance, half-width, and test results.
}
\examples{
# Generate sample MCMC chains
library(goftest)
chain1 <- rnorm(1000, mean = 1, sd = 4)
chain2 <- rnorm(1000, mean = 2, sd = 3)
mcmc_multi <- cbind(chain1, chain2)
mcmc_obj <- coda::mcmc(mcmc_multi)

# Run advheidel.diag with default settings (Kolmogorov-Smirnov test)
result_default <- advheidel.diag(mcmc_obj)
print(result_default)
summary(result_default)

# Run advheidel.diag with custom metrics enabled
result_multi <- advheidel.diag(mcmc_obj, custommetrics = TRUE)
print(result_multi)

# Use the Cramer-von Mises (CVM) test instead of Kolmogorov-Smirnov (KS)
result_cvm <- advheidel.diag(mcmc_obj, test_type = "cvm", custommetrics = TRUE)
print(result_cvm)

# Use the Shapiro-Wilk test
result_shapiro <- advheidel.diag(mcmc_obj, test_type = "shapiro", custommetrics = TRUE)
print(result_shapiro)

# Plot the results for one of the tests
plot(result_cvm)

}
