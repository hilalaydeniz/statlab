#' Frequency table of a variable (ni, fi, Ni, Fi)
#' @param x Vector of observations (numeric, character or factor).
#' @return Data frame with absolute (ni), relative (fi) and cumulative (Ni, Fi) frequencies.
#' @examples
#' freq_table(c(6, 8, 6, 8, 12, 6))
#' @export
freq_table <- function(x) {
  tab <- table(x)
  ni <- as.numeric(tab)
  N <- sum(ni)
  vals <- names(tab)
  if (is.numeric(x)) vals <- as.numeric(vals)
  data.frame(x = vals, ni = ni, fi = ni / N, Ni = cumsum(ni),
             Fi = cumsum(ni) / N, stringsAsFactors = FALSE)
}

#' Frequency table with class intervals (L_{i-1}, L_i]
#' @param x Numeric vector.
#' @param breaks Interval limits, e.g. c(0, 10, 20, 30).
#' @return Data frame with class mark, width, frequencies and density.
#' @examples
#' grouped_freq_table(c(1, 4, 7, 12, 15, 22), breaks = c(0, 10, 20, 30))
#' @export
grouped_freq_table <- function(x, breaks) {
  cls <- cut(x, breaks = breaks, include.lowest = TRUE)
  ni <- as.numeric(table(cls))
  N <- sum(ni)
  L0 <- head(breaks, -1)
  L1 <- breaks[-1]
  data.frame(interval = levels(cls), lower = L0, upper = L1,
             mark = (L0 + L1) / 2, width = L1 - L0, ni = ni, fi = ni / N,
             Ni = cumsum(ni), Fi = cumsum(ni) / N,
             density = ni / (L1 - L0), stringsAsFactors = FALSE)
}

#' Mean and variance from grouped data (class marks and frequencies)
#' @param mark Class marks.
#' @param ni Absolute frequencies.
#' @return List with mean, variance (population) and sd.
#' @export
grouped_stats <- function(mark, ni) {
  N <- sum(ni)
  m <- sum(mark * ni) / N
  v <- sum((mark - m)^2 * ni) / N
  list(mean = m, variance = v, sd = sqrt(v))
}

#' Mode(s) of a variable
#' @param x Vector.
#' @return All values with maximal frequency, or NA if there is no mode.
#' @export
mode_stat <- function(x) {
  tab <- table(x)
  m <- names(tab)[tab == max(tab)]
  if (length(tab) > 1 && length(m) == length(tab)) return(NA)
  if (is.numeric(x)) as.numeric(m) else m
}

#' Population variance (divides by N)
#' @param x Numeric vector.
#' @export
variance_pop <- function(x) mean((x - mean(x))^2)

#' Population standard deviation
#' @param x Numeric vector.
#' @export
sd_pop <- function(x) sqrt(variance_pop(x))

#' Sample (quasi-)variance (divides by N-1)
#' @param x Numeric vector.
#' @export
variance_sample <- function(x) sum((x - mean(x))^2) / (length(x) - 1)

#' Pearson coefficient of variation
#' @param x Numeric vector.
#' @export
cv <- function(x) sd_pop(x) / abs(mean(x))

#' Central moment of order r
#' @param x Numeric vector.
#' @param r Order.
#' @export
central_moment <- function(x, r) mean((x - mean(x))^r)

#' Fisher skewness coefficient g1 = m3 / s^3
#' @param x Numeric vector.
#' @export
skewness_fisher <- function(x) central_moment(x, 3) / sd_pop(x)^3

#' Pearson skewness coefficient (mean - mode) / s
#' @param x Numeric vector.
#' @export
skewness_pearson <- function(x) (mean(x) - mode_stat(x)[1]) / sd_pop(x)

#' Kurtosis coefficient g2 = m4 / s^4 - 3
#' @param x Numeric vector.
#' @export
kurtosis_fisher <- function(x) central_moment(x, 4) / sd_pop(x)^4 - 3

#' Quantiles using the course convention
#' @description Position N*p. If it is an integer the value is the average of
#' the observations at that position and the next one; otherwise the next
#' integer position is used.
#' @param x Numeric vector.
#' @param p Probabilities in (0, 1), e.g. c(0.25, 0.5, 0.75).
#' @examples
#' quantile_course(c(1, 3, 5, 7, 9, 11), c(0.25, 0.5, 0.75))
#' @export
quantile_course <- function(x, p) {
  x <- sort(x)
  N <- length(x)
  vapply(p, function(pp) {
    pos <- N * pp
    if (abs(pos - round(pos)) < 1e-9) {
      pos <- round(pos)
      if (pos < 1) x[1] else if (pos >= N) x[N] else (x[pos] + x[pos + 1]) / 2
    } else {
      x[ceiling(pos)]
    }
  }, numeric(1))
}

#' Outliers by the interquartile-range rule
#' @param x Numeric vector.
#' @param k Fence multiplier (1.5 mild, 3 extreme).
#' @return List with quartiles, fences and the outlying values.
#' @export
outliers_iqr <- function(x, k = 1.5) {
  q <- quantile_course(x, c(0.25, 0.75))
  iqr <- q[2] - q[1]
  lo <- q[1] - k * iqr
  hi <- q[2] + k * iqr
  list(Q1 = q[1], Q3 = q[2], IQR = iqr, lower_fence = lo, upper_fence = hi,
       outliers = x[x < lo | x > hi])
}

#' Standardise (z-scores) using the population sd
#' @param x Numeric vector.
#' @export
zscore <- function(x) (x - mean(x)) / sd_pop(x)

#' Rescale to [0, 1]
#' @param x Numeric vector.
#' @export
min_max <- function(x) (x - min(x)) / (max(x) - min(x))

#' One-shot summary of a numeric variable
#' @param x Numeric vector.
#' @return Named numeric vector with the main descriptive measures.
#' @examples
#' describe(c(2, 4, 4, 5, 9))
#' @export
describe <- function(x) {
  q <- quantile_course(x, c(0.25, 0.5, 0.75))
  c(n = length(x), mean = mean(x), median = q[2], Q1 = q[1], Q3 = q[3],
    min = min(x), max = max(x), range = max(x) - min(x),
    variance = variance_pop(x), sd = sd_pop(x), cv = cv(x),
    skewness = skewness_fisher(x), kurtosis = kurtosis_fisher(x))
}
