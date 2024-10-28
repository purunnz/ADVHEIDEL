#include <Rcpp.h>
#include <numeric>  // for accumulate
#include <cmath>    // for sqrt

using namespace Rcpp;

//' Heidelberger-Welch Metric Calculation with Stepwise Discard Ratios and Stationarity Test
 //'
 //' This function calculates the mean, variance, and half-width for a confidence interval
 //' on a given numeric vector at different discard levels. For each discard level (0%, 10%, ... 50%),
 //' it performs the Stationarity Test (Kolmogorov-Smirnov, Cramer-von Mises, or Shapiro-Wilk)
 //' and Half-Width Test, and returns the results.
 //'
 //' @param x A numeric vector containing the sample data.
 //' @param alpha A double representing the significance level for the confidence interval (default is 0.05).
 //' @param eps A double threshold for the half-width ratio test (default is 0.1).
 //' @param test_type Character. The stationarity test to apply. Options are "ks" for Kolmogorov-Smirnov,
 //'                  "cvm" for Cramer-von Mises, "shapiro" for Shapiro-Wilk test. Default is "ks".
 //' @return A DataFrame containing the results for each discard level, including mean, variance, p-value,
 //'         half-width ratio, test results, and the start iteration.
 //' @examples
 //' library(goftest)
 //' x <- rnorm(1000, mean = 0, sd = 1)
 //' result <- heidelWelchMetric(x)
 //' print(result)
 //'
 //' @export
 // [[Rcpp::export]]
 DataFrame heidelWelchMetric(NumericVector x, double alpha = 0.05, double eps = 0.1, std::string test_type = "ks") {
   int n = x.size();
   if (n < 1) {
     stop("Input vector must contain at least one value.");
   }

   // Check if `goftest` is loaded
   Rcpp::Function requireNamespace("requireNamespace");
   if (!Rcpp::as<bool>(requireNamespace("goftest", Named("quietly") = true))) {
     stop("Package 'goftest' is required but not installed.");
   }

   // Define the discard ratios (0%, 10%, ... 50%) and initialize the vectors to store results
   NumericVector discard_ratios = NumericVector::create(0, 0.1, 0.2, 0.3, 0.4, 0.5);
   NumericVector means, variances, halfwidths, halfwidth_ratios, p_values, start_iterations;
   LogicalVector halfwidth_passed, stationarity_passed;

   // Define R functions for the tests
   Rcpp::Function cvm_test("cvm.test"), ks_test("ks.test"), shapiro_test("shapiro.test");

   for (double discard_ratio : discard_ratios) {
     // Calculate number of elements to discard and the start iteration
     int discard_n = static_cast<int>(n * discard_ratio);
     int start_iteration = discard_n + 1;
     start_iterations.push_back(start_iteration);

     // Trim the vector based on the discard level
     NumericVector x_trimmed = x[Range(discard_n, n - 1)];
     int n_trimmed = x_trimmed.size();

     double mean = std::accumulate(x_trimmed.begin(), x_trimmed.end(), 0.0) / n_trimmed;
     means.push_back(mean);

     if (n_trimmed == 1) {
       variances.push_back(NA_REAL);
       halfwidths.push_back(NA_REAL);
       halfwidth_ratios.push_back(NA_REAL);
       halfwidth_passed.push_back(NA_LOGICAL);
       stationarity_passed.push_back(NA_LOGICAL);
       p_values.push_back(NA_REAL);
       continue;
     }

     double variance = 0.0;
     for (int i = 0; i < n_trimmed; ++i) {
       variance += std::pow(x_trimmed[i] - mean, 2);
     }
     variance /= (n_trimmed - 1);
     variances.push_back(variance);

     double t_value = R::qt(1 - alpha / 2, n_trimmed - 1, 1, 0);
     double halfwidth = t_value * std::sqrt(variance / n_trimmed);
     double halfwidth_ratio = std::abs(halfwidth / mean);
     bool hw_passed = halfwidth_ratio < eps;

     halfwidths.push_back(halfwidth);
     halfwidth_ratios.push_back(halfwidth_ratio);
     halfwidth_passed.push_back(hw_passed);

     // Perform stationarity test based on the test_type
     List test_result;
     double p_value;
     bool stationary;

     if (test_type == "cvm") {
       test_result = cvm_test(x_trimmed, Named("null") = "pnorm",
                              Named("mean") = mean, Named("sd") = std::sqrt(variance));
       p_value = as<double>(test_result["p.value"]);
     } else if (test_type == "shapiro") {
       test_result = shapiro_test(x_trimmed);
       p_value = as<double>(test_result["p.value"]);
     } else {  // Default to "ks" (Kolmogorov-Smirnov test)
       test_result = ks_test(x_trimmed, "pnorm", mean, std::sqrt(variance));
       p_value = as<double>(test_result["p.value"]);
     }

     stationary = p_value >= alpha;
     p_values.push_back(p_value);
     stationarity_passed.push_back(stationary);
   }

   // Return results in a DataFrame for easy viewing
   return DataFrame::create(
     Named("Discard Ratio") = discard_ratios,
     Named("Start Iteration") = start_iterations,
     Named("Stationarity Test") = stationarity_passed,
     Named("p-value") = p_values,
     Named("Mean") = means,
     Named("Variance") = variances,
     Named("Half-Width") = halfwidths,
     Named("Half-Width Ratio") = halfwidth_ratios,
     Named("Half-Width Test") = halfwidth_passed
   );
 }
