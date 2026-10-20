#' Moving average of order k
#' @param x Numeric series.
#' @param k Order.
#' @param centered For even k, apply the 2 x k centring so values align with t.
#' @return Numeric vector of length(x) with NA at the ends.
#' @examples
#' moving_average(c(3, 5, 7, 9, 11, 13), 3)
#' @export
moving_average <- function(x, k, centered = TRUE) {
  if (k %% 2 == 1 || !centered) {
    w <- rep(1 / k, k)
  } else {
    w <- c(0.5, rep(1, k - 1), 0.5) / k
  }
  as.numeric(stats::filter(x, w, sides = 2))
}

#' Least squares trend (polynomial in t)
#' @param x Numeric series.
#' @param t Time index (default 1..n).
#' @param degree 1 = line, 2 = parabola.
#' @return List with coefficients, fitted trend and a predict function.
#' @export
trend_ls <- function(x, t = seq_along(x), degree = 1) {
  d <- data.frame(x = x, t = t)
  m <- if (degree == 1) lm(x ~ t, d) else lm(x ~ t + I(t^2), d)
  list(coefficients = coef(m), fitted = as.numeric(fitted(m)),
       predict = function(tnew) as.numeric(predict(m, data.frame(t = tnew))))
}

#' Seasonal variation indices
#' @param x Numeric series.
#' @param period Number of seasons per cycle (e.g. 4 quarters, 12 months).
#' @param model "multiplicative" (ratios, mean 1) or "additive" (differences, mean 0).
#' @param trend "ma" (centred moving average) or "ls" (least squares line).
#' @return List with raw and corrected indices.
#' @export
seasonal_index <- function(x, period, model = c("multiplicative", "additive"),
                           trend = c("ma", "ls")) {
  model <- match.arg(model)
  trend <- match.arg(trend)
  tt <- seq_along(x)
  tr <- if (trend == "ma") moving_average(x, period, centered = TRUE) else trend_ls(x)$fitted
  raw <- if (model == "multiplicative") x / tr else x - tr
  season <- ((tt - 1) %% period) + 1
  idx <- as.numeric(tapply(raw, season, mean, na.rm = TRUE))
  corrected <- if (model == "multiplicative") idx / mean(idx) else idx - mean(idx)
  list(model = model, trend = tr, raw_index = idx, corrected_index = corrected)
}

#' Remove seasonality from a series
#' @param x Numeric series.
#' @param index Corrected seasonal indices (length = period).
#' @param model "multiplicative" or "additive".
#' @export
deseasonalize <- function(x, index, model = c("multiplicative", "additive")) {
  model <- match.arg(model)
  s <- index[((seq_along(x) - 1) %% length(index)) + 1]
  if (model == "multiplicative") x / s else x - s
}

#' Forecast h steps ahead: least squares trend x seasonal index
#' @param x Numeric series.
#' @param period Seasons per cycle.
#' @param h Number of future periods.
#' @param model "multiplicative" or "additive".
#' @param trend Method used to estimate indices: "ma" or "ls".
#' @return Numeric vector of forecasts.
#' @export
forecast_seasonal <- function(x, period, h = period, model = c("multiplicative", "additive"),
                              trend = c("ma", "ls")) {
  model <- match.arg(model)
  trend <- match.arg(trend)
  si <- seasonal_index(x, period, model, trend)$corrected_index
  des <- deseasonalize(x, si, model)
  n <- length(x)
  d <- data.frame(des = des, tt = seq_len(n))
  m <- lm(des ~ tt, d)
  tf <- (n + 1):(n + h)
  tr_f <- as.numeric(predict(m, data.frame(tt = tf)))
  s_f <- si[((tf - 1) %% period) + 1]
  if (model == "multiplicative") tr_f * s_f else tr_f + s_f
}

#' Simple index numbers I_t/t0 = 100 * x_t / x_t0
#' @param x Numeric series.
#' @param base Position of the base period.
#' @export
simple_index <- function(x, base = 1) 100 * x / x[base]

#' Change the base of an index series
#' @param index Index series.
#' @param new_base Position of the new base period.
#' @export
change_base <- function(index, new_base = 1) 100 * index / index[new_base]

#' Percentage variation between consecutive periods
#' @param x Numeric series.
#' @export
pct_change <- function(x) c(NA, diff(x) / head(x, -1) * 100)

#' Deflate a nominal series (real values in base-year units)
#' @param nominal Nominal values.
#' @param index Price index (base 100).
#' @export
deflate <- function(nominal, index) nominal / index * 100

#' Laspeyres price index
#' @param p0 Base prices.
#' @param pt Current prices.
#' @param q0 Base quantities.
#' @export
laspeyres <- function(p0, pt, q0) 100 * sum(pt * q0) / sum(p0 * q0)

#' Paasche price index
#' @param p0 Base prices.
#' @param pt Current prices.
#' @param qt Current quantities.
#' @export
paasche <- function(p0, pt, qt) 100 * sum(pt * qt) / sum(p0 * qt)

#' Fisher ideal index (geometric mean of Laspeyres and Paasche)
#' @param p0,pt,q0,qt Prices and quantities in base (0) and current (t) periods.
#' @export
fisher_index <- function(p0, pt, q0, qt) sqrt(laspeyres(p0, pt, q0) * paasche(p0, pt, qt))
