
# ADVHEIDEL

The **Heidelberger-Welch Test** is used to verify that MCMC samples have stabilized and are drawn from a consistent distribution, which is essential for reliable results, particularly in Bayesian analysis.

### Stationarity Test
The stationarity test checks if the MCMC chain has reached a stable state, ensuring the samples reflect a consistent distribution. This stability is critical for analyzing dependable data.

1. **Start with Full Data**: The test begins by examining the entire chain.
2. **Incremental Discards**: If instability is detected, the first 10% of samples are removed, and the test is repeated.
3. **Repeat Up to 50%**: The process continues, removing 10% increments (up to 50%) until stability is achieved.
4. **Use Stable Data**: Once stability is confirmed, only the stable part of the chain is used for further analysis.

This approach discards early, potentially unreliable samples, focusing only on the stable portion.

### Halfwidth Test
The halfwidth test ensures that the number of samples is sufficient for a reliable mean estimate by calculating a confidence interval around the mean to assess precision.

1. **Calculate Confidence Interval**: Using stable samples, it calculates the range where the true mean likely falls.
2. **Check Precision**: The width of this interval (halfwidth) is compared to the mean. A smaller halfwidth relative to the mean indicates precision.
3. **Pass or Fail**:
   - **Pass**: A small halfwidth confirms enough samples for a reliable mean.
   - **Fail**: A large halfwidth suggests that more samples may be needed for a precise result.

In summary:
- The **Stationarity Test** confirms that we’re working with stable data.
- The **Halfwidth Test** ensures that our mean estimate is accurate and dependable. 



>>>>>>> 2abc78c (first commit)
