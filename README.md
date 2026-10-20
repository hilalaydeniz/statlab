# Statistics course repository

| Folder | File | Contents |
|---|---|---|
| 01_descriptive_statistics/R | descriptive.R | frequency tables, mean/mode/quantiles, variance, skewness, kurtosis, outliers |
| 02_statistical_models/R | models.R | joint/marginal/conditional tables, covariance, correlation, regression, transformations |
| 03_time_series/R | timeseries.R | moving averages, trend, seasonal indices, forecast, index numbers |
| 04_probability/R | probability.R | classical probability, conditional, total probability, Bayes, diagnostic tests |
| 05_random_variables/R | distributions.R | discrete uniform, Poisson, binomial, geometric, hypergeometric (pmf/cdf/rnd/moments) |

## Use
```r
source("load_all.R")
bayes(c(0.2, 0.8), c(0.9, 0.4))
```
Note: models.R uses `sd_pop()` and `variance_pop()` from descriptive.R, so load everything with `load_all.R`.
