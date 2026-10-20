#' Discrete Uniform Distribution U(a,b)
#' @family discrete-distributions
#' @description PMF, CDF, RNG and moments for a discrete uniform distribution on integers a..b.
#' @param k Value(s) at which to evaluate the PMF.
#' @param q Thresholds at which to evaluate the CDF.
#' @param a Integer lower bound.
#' @param b Integer upper bound (b >= a).
#' @param n Non-negative integer sample size.
#' @examples
#' discrete_uniform_pmf(1:6, 1, 6)
#' discrete_uniform_cdf(c(0, 3, 6), 1, 6)
#' discrete_uniform_rnd(5, 1, 6)
#' discrete_uniform_moments(1, 6)
#' @export
discrete_uniform_pmf <- function(k, a, b) {
  if (any(!is.finite(c(a, b))) || a %% 1 != 0 || b %% 1 != 0 || b < a) {
    stop("'a' and 'b' must be finite integers with b >= a.")
  }
  k <- as.numeric(k)
  ok <- (k %% 1 == 0) & (k >= a) & (k <= b)
  out <- numeric(length(k))
  out[ok] <- 1 / (b - a + 1)
  out
}

#' @rdname discrete_uniform_pmf
#' @export
discrete_uniform_cdf <- function(q, a, b) {
  if (any(!is.finite(c(a, b))) || a %% 1 != 0 || b %% 1 != 0 || b < a) {
    stop("'a' and 'b' must be finite integers with b >= a.")
  }
  q <- as.numeric(q)
  Fq <- numeric(length(q))
  Fq[q >= b] <- 1
  inside <- (q >= a) & (q < b)
  Fq[inside] <- (floor(q[inside]) - a + 1) / (b - a + 1)
  Fq
}

#' @rdname discrete_uniform_pmf
#' @export
discrete_uniform_rnd <- function(n, a, b) {
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 0 || n %% 1 != 0) {
    stop("'n' must be a non-negative integer.")
  }
  if (any(!is.finite(c(a, b))) || a %% 1 != 0 || b %% 1 != 0 || b < a) {
    stop("'a' and 'b' must be finite integers with b >= a.")
  }
  as.integer(a - 1 + sample.int(b - a + 1, size = n, replace = TRUE))
}

#' @rdname discrete_uniform_pmf
#' @export
discrete_uniform_moments <- function(a, b) {
  if (any(!is.finite(c(a, b))) || a %% 1 != 0 || b %% 1 != 0 || b < a) {
    stop("'a' and 'b' must be finite integers with b >= a.")
  }
  list(mean = (a + b) / 2, variance = ((b - a + 1)^2 - 1) / 12)
}

#' Poisson Distribution P(lambda)
#' @family discrete-distributions
#' @description PMF, CDF, RNG and moments of a Poisson distribution, lambda > 0.
#' @param k Non-negative integer value(s).
#' @param q Thresholds for the CDF.
#' @param lambda Positive rate.
#' @param n Sample size.
#' @examples
#' poisson_pmf(0:5, 3)
#' poisson_cdf(0:5, 3)
#' poisson_moments(3)
#' @export
poisson_pmf <- function(k, lambda) {
  if (!is.numeric(lambda) || any(lambda <= 0)) stop("'lambda' must be positive.")
  k <- as.numeric(k)
  if (any(k < 0 | k %% 1 != 0)) stop("'k' must be non-negative integers.")
  dpois(k, lambda)
}

#' @rdname poisson_pmf
#' @export
poisson_cdf <- function(q, lambda) {
  if (!is.numeric(lambda) || any(lambda <= 0)) stop("'lambda' must be positive.")
  ppois(as.numeric(q), lambda)
}

#' @rdname poisson_pmf
#' @export
poisson_rnd <- function(n, lambda) {
  if (!is.numeric(lambda) || any(lambda <= 0)) stop("'lambda' must be positive.")
  if (!is.numeric(n) || length(n) != 1 || n < 0 || n %% 1 != 0) stop("'n' must be a non-negative integer.")
  rpois(n, lambda)
}

#' @rdname poisson_pmf
#' @export
poisson_moments <- function(lambda) {
  if (!is.numeric(lambda) || any(lambda <= 0)) stop("'lambda' must be positive.")
  list(mean = lambda, variance = lambda)
}

#' Binomial Distribution B(size, prob)
#' @family discrete-distributions
#' @description PMF, CDF, RNG and moments of a Binomial distribution.
#' @param k Number of successes.
#' @param q Thresholds for the CDF.
#' @param size Number of trials (non-negative integer).
#' @param prob Success probability in [0, 1].
#' @param n Sample size.
#' @examples
#' binomial_pmf(0:4, 4, 0.5)
#' binomial_cdf(2, 4, 0.5)
#' binomial_moments(10, 0.3)
#' @export
binomial_pmf <- function(k, size, prob) {
  .check_binom(size, prob)
  dbinom(as.numeric(k), size, prob)
}

#' @rdname binomial_pmf
#' @export
binomial_cdf <- function(q, size, prob) {
  .check_binom(size, prob)
  pbinom(as.numeric(q), size, prob)
}

#' @rdname binomial_pmf
#' @export
binomial_rnd <- function(n, size, prob) {
  .check_binom(size, prob)
  rbinom(n, size, prob)
}

#' @rdname binomial_pmf
#' @export
binomial_moments <- function(size, prob) {
  .check_binom(size, prob)
  list(mean = size * prob, variance = size * prob * (1 - prob))
}

.check_binom <- function(size, prob) {
  if (!is.numeric(size) || size < 0 || size %% 1 != 0) stop("'size' must be a non-negative integer.")
  if (!is.numeric(prob) || prob < 0 || prob > 1) stop("'prob' must be in [0, 1].")
}

#' Geometric Distribution (number of trials until the first success)
#' @family discrete-distributions
#' @description Support {1, 2, 3, ...}: P(X = k) = (1 - p)^(k - 1) p.
#' Note: R's dgeom() counts failures, so k - 1 is used internally.
#' @param k Trial number of the first success (>= 1).
#' @param q Thresholds for the CDF.
#' @param prob Success probability in (0, 1].
#' @param n Sample size.
#' @examples
#' geometric_pmf(1:5, 0.3)
#' geometric_moments(0.3)
#' @export
geometric_pmf <- function(k, prob) {
  .check_geom(prob)
  k <- as.numeric(k)
  out <- numeric(length(k))
  ok <- k >= 1 & k %% 1 == 0
  out[ok] <- (1 - prob)^(k[ok] - 1) * prob
  out
}

#' @rdname geometric_pmf
#' @export
geometric_cdf <- function(q, prob) {
  .check_geom(prob)
  q <- floor(as.numeric(q))
  ifelse(q < 1, 0, 1 - (1 - prob)^q)
}

#' @rdname geometric_pmf
#' @export
geometric_rnd <- function(n, prob) {
  .check_geom(prob)
  rgeom(n, prob) + 1
}

#' @rdname geometric_pmf
#' @export
geometric_moments <- function(prob) {
  .check_geom(prob)
  list(mean = 1 / prob, variance = (1 - prob) / prob^2)
}

.check_geom <- function(prob) {
  if (!is.numeric(prob) || prob <= 0 || prob > 1) stop("'prob' must be in (0, 1].")
}

#' Hypergeometric Distribution (sampling without replacement)
#' @family discrete-distributions
#' @description Population of N elements with K successes; a sample of size n is drawn.
#' @param k Number of successes in the sample.
#' @param q Thresholds for the CDF.
#' @param N Population size.
#' @param K Number of successes in the population.
#' @param n Sample size.
#' @examples
#' # 20 components, 3 defective, inspect 4: P(exactly 1 defective)
#' hypergeometric_pmf(1, N = 20, K = 3, n = 4)  # 0.42
#' @export
hypergeometric_pmf <- function(k, N, K, n) dhyper(k, K, N - K, n)

#' @rdname hypergeometric_pmf
#' @export
hypergeometric_cdf <- function(q, N, K, n) phyper(q, K, N - K, n)

#' @rdname hypergeometric_pmf
#' @export
hypergeometric_moments <- function(N, K, n) {
  p <- K / N
  list(mean = n * p, variance = n * p * (1 - p) * (N - n) / (N - 1))
}

