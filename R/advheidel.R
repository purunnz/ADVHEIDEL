#' advheidel.diag Function: Enhanced Heidelberger-Welch MCMC Diagnostics
#'
#' This function performs an advanced Heidelberger-Welch MCMC diagnostic on one or more chains
#' with options for different stationarity tests. It includes the Kolmogorov-Smirnov (KS), Cramer-von Mises (CVM),
#' Shapiro-Wilk tests, and custom metrics when applicable. The function also calculates the Half-Width ratio for each chain.
#'
#' @param mcmc_chain An MCMC object (matrix or data frame) containing one or more chains.
#' @param alpha Numeric. The significance level for the stationarity and half-width tests. Default is 0.05.
#' @param eps Numeric. The threshold for the Half-Width ratio test. Default is 0.1.
#' @param test_type Character. The stationarity test to apply. Options are "ks" for Kolmogorov-Smirnov,
#'        "cvm" for Cramer-von Mises, "shapiro" for Shapiro-Wilk test. Default is "ks".
#' @param custommetrics Logical. If TRUE, the function calculates additional custom metrics using the
#'        \code{heidelWelchMetric} function. Default is FALSE.
#'
#' @details
#' This function applies the Heidelberger-Welch stationarity and Half-Width tests to diagnose convergence issues in MCMC chains.
#' It allows users to choose from three different stationarity tests (\code{ks}, \code{cvm}, \code{shapiro}) and provides additional metrics such as the Half-Width ratio.
#' The \code{custommetrics} parameter, if set to TRUE, enables the calculation of additional diagnostic metrics using the C++ implementation via the \code{heidelWelchMetric} function.
#'
#' The stationarity test involves discarding parts of the MCMC chains and checking for convergence. If custom metrics are enabled, the function returns diagnostic results, including mean, variance, half-width, and test results.
#'
#' @importFrom graphics abline legend par plot text
#' @importFrom stats sd qt ks.test shapiro.test
#' @importFrom goftest cvm.test
#'
#' @return
#' A list object of class \code{advHeidelDiag} containing:
#' \describe{
#'   \item{results}{A data frame with diagnostic results for each chain.}
#'   \item{mcmc_chain}{The input MCMC chain.}
#' }
#'
#' @examples
#' # Generate sample MCMC chains
#' library(goftest)
#' chain1 <- rnorm(1000, mean = 1, sd = 4)
#' chain2 <- rnorm(1000, mean = 2, sd = 3)
#' mcmc_multi <- cbind(chain1, chain2)
#' mcmc_obj <- coda::mcmc(mcmc_multi)
#'
#' # Run advheidel.diag with default settings (Kolmogorov-Smirnov test)
#' result_default <- advheidel.diag(mcmc_obj)
#' print(result_default)
#' summary(result_default)
#'
#' # Run advheidel.diag with custom metrics enabled
#' result_multi <- advheidel.diag(mcmc_obj, custommetrics = TRUE)
#' print(result_multi)
#'
#' # Use the Cramer-von Mises (CVM) test instead of Kolmogorov-Smirnov (KS)
#' result_cvm <- advheidel.diag(mcmc_obj, test_type = "cvm", custommetrics = TRUE)
#' print(result_cvm)
#'
#' # Use the Shapiro-Wilk test
#' result_shapiro <- advheidel.diag(mcmc_obj, test_type = "shapiro", custommetrics = TRUE)
#' print(result_shapiro)
#'
#' # Plot the results for one of the tests
#' plot(result_cvm)
#'
#' @export
advheidel.diag <- function(mcmc_chain, alpha = 0.05, eps = 0.1, test_type = "ks", custommetrics = FALSE) {

  # ---------------------------------------------------------------
  # Early check: Reject list inputs that are not data frames
  if (is.list(mcmc_chain) && !is.data.frame(mcmc_chain)) {
    stop("Input 'mcmc_chain' must be a matrix, data frame, or MCMC object.")
  }

  # Validate test_type input
  valid_test_types <- c("ks", "cvm", "shapiro")
  if (!test_type %in% valid_test_types) {
    stop("Unsupported test type. Choose 'ks', 'cvm', or 'shapiro'.")
  }

  # Convert single vector input to a matrix with one column immediately
  if (is.vector(mcmc_chain)) {
    mcmc_chain <- matrix(mcmc_chain, ncol = 1)
  }


  # Early check: Ensure mcmc_chain is either a matrix, data frame, or MCMC object
  if (!is.matrix(mcmc_chain) && !inherits(mcmc_chain, "mcmc") && !is.data.frame(mcmc_chain)) {
    stop("Input 'mcmc_chain' must be a matrix, data frame, or MCMC object.")
  }

  # Check if input is an invalid list
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("Error: 'alpha' must be a numeric value between 0 and 1.")
  }
  if (!is.numeric(eps) || eps <= 0) {
    stop("Error: 'eps' must be a positive numeric value.")
  }

  # Convert to MCMC object if not already
  if (!inherits(mcmc_chain, "mcmc")) {
    mcmc_chain <- coda::as.mcmc(mcmc_chain)
  }

  # Ensure mcmc_chain has columns
  if (is.null(ncol(mcmc_chain))) {
    mcmc_chain <- matrix(mcmc_chain, ncol = 1)
  }

  # Convert single vector input to a matrix with one column before type validation
  if (is.vector(mcmc_chain)) {
    mcmc_chain <- matrix(mcmc_chain, ncol = 1)
  }



  # Check for null input
  if (is.null(mcmc_chain)) {
    stop("Error: The 'mcmc_chain' is NULL. Please provide valid data.")
  }

  # Check if mcmc_chain is a list and convert if needed
  if (is.list(mcmc_chain)) {
    if (all(sapply(mcmc_chain, is.numeric))) {
      # Convert list of numeric vectors to matrix
      mcmc_chain <- do.call(cbind, mcmc_chain)
    } else {
      stop("Input 'mcmc_chain' must be a matrix, data frame, or MCMC object.")
    }
  }

  # Convert data frame to matrix if necessary
  if (is.data.frame(mcmc_chain)) {
    mcmc_chain <- as.matrix(mcmc_chain)
  }

  # Convert single vector input to a matrix with one column
  if (is.vector(mcmc_chain)) {
    mcmc_chain <- matrix(mcmc_chain, ncol = 1)
  }

  # Validate the input type after conversion
  if (!is.matrix(mcmc_chain) && !inherits(mcmc_chain, "mcmc")) {
    stop("Input 'mcmc_chain' must be a matrix, data frame, or MCMC object.")
  }

  # Check for invalid values (Inf, NA, NaN)
  if (any(is.infinite(mcmc_chain) | is.na(mcmc_chain) | is.nan(mcmc_chain))) {
    stop("Error: The 'mcmc_chain' contains Inf, NA, or NaN values. Please provide valid data.")
  }

  # Convert to MCMC object if not already
  if (!inherits(mcmc_chain, "mcmc")) {
    mcmc_chain <- coda::as.mcmc(mcmc_chain)
  }

  # Apply diagnostic tests to each chain
  results <- lapply(1:ncol(mcmc_chain), function(chain_index) {
    chain <- as.numeric(mcmc_chain[, chain_index])

    # Run Heidelberger-Welch metric, which includes both stationarity and half-width tests
    metric_result <- heidelWelchMetric(chain, alpha, eps, test_type)

    # # Debug output for metric_result
    # print(paste("Debug: Metric result for Chain", chain_index))
    # print(metric_result)

    # Find the first row where Stationarity Test is TRUE
    first_stationary_row_index <- which(as.logical(metric_result[["Stationarity.Test"]]))[1]

    # Debug output for first_stationary_row_index
    # print(paste("Debug: First stationary row index for Chain", chain_index, "is", first_stationary_row_index))

    # Check if any row passed the stationarity test
    if (!is.na(first_stationary_row_index)) {
      # Extract the first row that passes the stationarity test, keeping it as a data frame
      final_row <- metric_result[first_stationary_row_index, , drop = FALSE]
    } else {
      # Set default values if no stationarity test passed, matching the structure of metric_result
      # Extract the last row of metric_result if no stationarity test passed
      final_row <- metric_result[nrow(metric_result), , drop = FALSE]
      #
      # final_row <- data.frame(
      #   "Stationarity.Test" = FALSE,
      #   "Discard.Ratio" = NA,
      #   "p.value" = NA,
      #   "Mean" = NA,
      #   "Variance" = NA,
      #   "Half.Width" = NA,
      #   "Half.Width.Ratio" = NA,
      #   "Half.Width.Test" = NA
      # )
    }

    # Debug: Print column names and data of final_row
    # print("Debug: final_row content")
    # print(final_row)
    # print("Debug: final_row column names")
    # print(names(final_row))


    # Compile results using values from the first stationary row or default values
    result <- list(
      Chain = paste0("Chain", chain_index),
      StationarityTest = ifelse(final_row[["Stationarity.Test"]][1], "Passed", "Failed"),
      # StartIteration = if (!is.na(first_stationary_row_index)) first_stationary_row_index else NA,
      StartIteration = final_row[["Start.Iteration"]][1],
      PValue = final_row[["p.value"]][1],
      TestType = test_type,
      HalfWidthTest = ifelse(final_row[["Half.Width.Test"]][1], "Passed", "Failed"),
      Mean = final_row[["Mean"]][1],
      HalfWidth = final_row[["Half.Width"]][1],
      HalfWidthRatio = final_row[["Half.Width.Ratio"]][1],
      CustomMetric = if (custommetrics) metric_result else NULL
    )

    # Debug: Print result list to verify contents
    # print("Debug: Compiled result list")
    # print(result)

    return(result)
  })

  # Convert only primary diagnostic results (without CustomMetric) to data frame
  main_results <- do.call(rbind, lapply(results, function(x) {
    as.data.frame(x[!names(x) %in% "CustomMetric"], stringsAsFactors = FALSE)
  }))

  # Embed CustomMetric lists for each chain to retain the structure
  main_results$CustomMetric <- lapply(results, function(x) x$CustomMetric)

  # Return the structured list object
  result_obj <- list(results = main_results, mcmc_chain = mcmc_chain)
  class(result_obj) <- "advHeidelDiag"
  return(result_obj)
}


# ---------------------------------------------------------------
#' Print Method for advHeidelDiag Class
#'
#' This function prints the Heidelberger-Welch diagnostic results in a table format,
#' including custom metrics and test type if available.
#'
#' @param x An object of class `advHeidelDiag` containing diagnostic results.
#' @param ... Additional arguments (not used).
#'
#' @method print advHeidelDiag
#' @export
print.advHeidelDiag <- function(x, ...) {
  # Check if the object contains the diagnostic results
  if (!"results" %in% names(x) || is.null(x$results)) {
    stop("The object does not contain the diagnostic results.")
  }

  diag_results <- x$results

  # Print header
  cat("Heidelberger-Welch MCMC Diagnostic Results:\n")
  cat(strrep("-", 140), "\n")

  # Check if custom metrics are present
  has_custom_metrics <- any(sapply(diag_results$CustomMetric, is.data.frame))

  # Adjust the table header based on whether custom metrics exist

  cat(sprintf("%-10s %-15s %-10s %-15s %-10s %-10s %-10s %-10s %-10s\n",
                "Chain", "Stationarity", "p-value", "Start Iteration", "TestType",
                "Mean", "Half-Width", "(Ratio)", "Half-Width Test"))

  cat(strrep("-", 140), "\n")

  # Loop through each chain and print results in one row per chain
  for (i in 1:nrow(diag_results)) {
    # Access main diagnostic results
    cat(sprintf("%-10d %-15s %-10.4f %-15d %-10s %-10.4f %-10.4f %-10.4f %-15s",
                i,
                diag_results$StationarityTest[i],
                diag_results$PValue[i],
                diag_results$StartIteration[i],
                diag_results$TestType[i],
                diag_results$Mean[i],
                diag_results$HalfWidth[i],
                diag_results$HalfWidthRatio[i],
                diag_results$HalfWidthTest[i]))
    cat("\n")
  }

  cat(strrep("-", 140), "\n")


  # Loop through each chain and print results in one row per chain
  for (i in 1:nrow(diag_results)) {

    # Check and print custom metrics if available for the chain
    if (has_custom_metrics && !is.null(diag_results$CustomMetric[[i]])) {
      custom_metrics <- diag_results$CustomMetric[[i]]

      # Print header for custom metrics
      cat(sprintf("\n  Custom Metrics for %s: \n", diag_results$Chain[i]))

      # Print header for custom metrics
      cat("  ", sprintf("%-15s %-15s %-10s %-15s %-10s %-10s %-15s %-10s %-10s\n",
                        "Discard Ratio", "Start.Iteration", "Stationary", "p-value", "Mean", "Variance", "Half-Width", "Ratio", "Test"))
      cat("  ", strrep("-", 120), "\n")

      # Print each row of custom metrics for the chain
      for (j in 1:nrow(custom_metrics)) {
        cat("  ", sprintf("%-15.2f %-15s %-10s %-15.4f %-10.4f %-10.4f %-15.4f %-10.4f %-10s\n",
                          custom_metrics$Discard.Ratio[j],
                          custom_metrics$Start.Iteration[j],
                          ifelse(custom_metrics$Stationarity.Test[j], "Passed", "Failed"),
                          custom_metrics$p.value[j],
                          custom_metrics$Mean[j],
                          custom_metrics$Variance[j],
                          custom_metrics$Half.Width[j],
                          custom_metrics$Half.Width.Ratio[j],
                          ifelse(custom_metrics$Half.Width.Test[j], "Passed", "Failed")
        ))
      }
      # Separator after each custom metric block
      cat(strrep("-", 140), "\n")
    } else {
      # If no custom metrics, just add a newline for spacing
      cat("\n")
    }
  }

  # Final separator line for the entire output
  cat(strrep("-", 140), "\n")


  # Print a summary at the end
  total_chains <- nrow(diag_results)
  passed_stationarity <- sum(diag_results$StationarityTest == "Passed", na.rm = TRUE)
  passed_halfwidth <- sum(diag_results$HalfWidthTest == "Passed", na.rm = TRUE)

  cat("\nSummary:\n")
  cat("  - Total Chains:                  ", total_chains, "\n")
  cat("  - Chains Passed Stationarity:    ", passed_stationarity, "/", total_chains, "\n")
  cat("  - Chains Passed Half-Width Test: ", passed_halfwidth, "/", total_chains, "\n")

  if (has_custom_metrics) {
    cat("  - Custom Metrics Available\n")
  }
}

# ---------------------------------------------------------------
#' Summary Method for advHeidelDiag Class
#'
#' This function provides a concise summary of the Heidelberger-Welch diagnostic results,
#' including the number of chains that passed or failed the stationarity and half-width tests.
#'
#' @param object An object of class `advHeidelDiag`.
#' @param ... Additional arguments (not used).
#'
#' @method summary advHeidelDiag
#' @export
summary.advHeidelDiag <- function(object, ...) {
  total_chains <- nrow(object$results)
  stationarity_passed <- sum(object$results$StationarityTest == "Passed", na.rm = TRUE)
  stationarity_failed <- sum(object$results$StationarityTest == "Failed", na.rm = TRUE)
  half_width_passed <- sum(object$results$HalfWidthTest == "Passed", na.rm = TRUE)
  half_width_failed <- sum(object$results$HalfWidthTest == "Failed", na.rm = TRUE)

  cat("Summary of Heidelberger-Welch Diagnostic Results:\n")
  cat("Total Chains Analyzed: ", total_chains, "\n")
  cat("Stationarity Test: \n")
  cat("  Passed: ", stationarity_passed, "\n")
  cat("  Failed: ", stationarity_failed, "\n")
  cat("Half-Width Test: \n")
  cat("  Passed: ", half_width_passed, "\n")
  cat("  Failed: ", half_width_failed, "\n")
}

#' Plot Method for advHeidelDiag Class
#'
#' This function generates diagnostic plots (trace, cumulative mean, half-width) for each chain.
#'
#' @param x An object of class `advHeidelDiag` containing MCMC data and diagnostic results.
#' @param ... Additional arguments (not used).
#'
#' @method plot advHeidelDiag
#' @export
plot.advHeidelDiag <- function(x, ...) {
  # Extract MCMC object and diagnostic results from the input object
  if (!"mcmc_chain" %in% names(x) || is.null(x$mcmc_chain)) {
    stop("The object does not contain the MCMC chain data.")
  }

  mcmc_obj <- x$mcmc_chain
  hw_diag_results <- x$results

  # Check if the MCMC object is correctly formatted
  if (!is.matrix(mcmc_obj)) {
    if (is.vector(mcmc_obj)) {
      mcmc_obj <- matrix(mcmc_obj, ncol = 1)  # Convert vector to matrix if necessary
    } else {
      stop("MCMC data must be either a vector or a matrix.")
    }
  }

  # Cumulative mean plot function
  cumulative_mean_plot <- function(chain_values, chain_index, chain_mean, p_value, start_iteration, convergence_label, test_type, stationarity_test_result) {
    iterations <- 1:length(chain_values)
    cumulative_means <- cumsum(chain_values) / iterations

    # Plot cumulative mean against iterations
    plot(iterations, cumulative_means, type = "l", col = "blue",
         xlab = "Iteration", ylab = "Cumulative Mean",
         main = paste("Cumulative Mean Plot for Chain", chain_index))

    # Add horizontal line for overall mean
    abline(h = chain_mean, col = "red", lty = 2)

    # Mark the start iteration with a vertical line
    abline(v = start_iteration, col = "green", lty = 2)
    text(x = start_iteration, y = max(cumulative_means, na.rm = TRUE),
         labels = paste("Start Iteration:", start_iteration), col = "green", pos = 4)

    # Update legend to show the test results
    legend("topright", legend = c("Cumulative Mean", "Overall Mean",
                                  paste("Stationarity Test:", stationarity_test_result),
                                  if (!is.na(p_value)) paste("p-value:", round(p_value, 3)) else "",
                                  paste("Start Iteration:", start_iteration)),
           col = c("blue", "red", "black", "black", "green"),
           lty = c(1, 2, NA, NA, 2), bg = "white", cex = 0.8)
  }

  # Half-width plot function
  plot_halfwidth <- function(chain_values, chain_index, mean_halfwidth, chain_mean, halfwidth_ratio, convergence_label, halfwidth_test_result) {
    iterations <- 1:length(chain_values)
    sd_values <- sapply(1:length(chain_values), function(i) {
      if (i > 1) sd(chain_values[1:i], na.rm = TRUE) else NA
    })

    # Calculate half-width for each iteration (using t-value for 95% CI)
    half_widths <- qt(0.975, pmax(iterations - 1, 1)) * sd_values / sqrt(pmax(iterations, 1))

    # Plot half-width against iterations
    plot(iterations, half_widths, type = "l", col = "purple",
         xlab = "Iteration", ylab = "Half-Width of CI",
         main = paste("Half-Width Plot for Chain", chain_index))

    # Add horizontal line for the mean half-width
    abline(h = mean_halfwidth, col = "red", lty = 2)

    # Mark the overall chain mean with a horizontal line
    abline(h = chain_mean, col = "blue", lty = 2)
    text(x = 0.8 * max(iterations), y = chain_mean,
         labels = paste("Mean:", round(chain_mean, 3)), col = "blue", pos = 4)

    # Update legend to show the half-width test result and ratio
    legend("topright", legend = c("Half-Width", paste("Mean Half-Width:", round(mean_halfwidth, 3)),
                                  paste("Overall Mean:", round(chain_mean, 3)),
                                  paste("HalfWidth Ratio:", round(halfwidth_ratio, 3)),
                                  paste("Halfwidth Test:", halfwidth_test_result)),
           col = c("purple", "red", "blue", "black", "black"),
           lty = c(1, 2, 2, NA, NA), bg = "white", cex = 0.8)
  }

  # Loop over each chain and generate plots
  for (i in 1:ncol(mcmc_obj)) {
    par(mfrow = c(1, 2))  # Create a 1x2 layout for the current chain

    # Extract chain values and diagnostics
    chain_values <- mcmc_obj[, i]
    chain_mean <- as.numeric(hw_diag_results[i, "Mean"])
    p_value <- as.numeric(hw_diag_results[i, "PValue"])
    mean_halfwidth <- as.numeric(hw_diag_results[i, "HalfWidth"])
    start_iteration <- as.numeric(hw_diag_results[i, "StartIteration"])
    halfwidth_ratio <- as.numeric(hw_diag_results[i, "HalfWidthRatio"])
    test_type <- as.character(hw_diag_results[i, "TestType"])

    # Determine convergence labels
    stationarity_test_result <- as.character(hw_diag_results[i, "StationarityTest"])
    halfwidth_test_result <- as.character(hw_diag_results[i, "HalfWidthTest"])
    convergence_label <- ifelse(stationarity_test_result == "Passed", "Converged", "Not Converged")

    # Trace plot for the current chain
    plot(chain_values, type = "l", col = ifelse(convergence_label == "Not Converged", "red", "blue"),
         main = paste("Trace Plot for Chain", i),
         xlab = "Iteration", ylab = "Value")

    # Cumulative mean plot for the current chain, including diagnostics
    cumulative_mean_plot(chain_values, i, chain_mean, p_value, start_iteration, convergence_label, test_type, stationarity_test_result)

    # Half-width plot for the current chain, including diagnostics
    plot_halfwidth(chain_values, i, mean_halfwidth, chain_mean, halfwidth_ratio, convergence_label, halfwidth_test_result)
  }
}


