#' Two-way frequency table with optional margins
#' @param x,y Vectors of the same length.
#' @param relative If TRUE, relative frequencies fij.
#' @param margins If TRUE, add marginal totals.
#' @export
joint_table <- function(x, y, relative = FALSE, margins = TRUE) {
  tab <- table(x, y)
  if (relative) tab <- prop.table(tab)
  if (margins) tab <- addmargins(tab)
  tab
}

#' Conditional distribution X|Y=y or Y|X=x
#' @param x,y Vectors of the same length.
#' @param given Value of the conditioning variable.
#' @param on "y" gives X | Y = given; "x" gives Y | X = given.
#' @export
cond_dist <- function(x, y, given, on = c("y", "x")) {
  on <- match.arg(on)
  tab <- table(x, y)
  v <- if (on == "y") tab[, as.character(given)] else tab[as.character(given), ]
  v / sum(v)
}

#' Expected frequencies under independence, chi-square, C and Cramer's V
#' @param x,y Categorical vectors.
#' @return List with expected table, chi2, contingency coefficient and Cramer's V.
#' @export
independence_check <- function(x, y) {
  obs <- table(x, y)
  n <- sum(obs)
  expd <- outer(rowSums(obs), colSums(obs)) / n
  chi2 <- sum((obs - expd)^2 / expd)
  m <- min(nrow(obs), ncol(obs)) - 1
  list(expected = expd, chi2 = chi2,
       contingency_C = sqrt(chi2 / (chi2 + n)),
       cramer_V = sqrt(chi2 / (n * m)))
}

#' Ordinary or central moment of order (r, s)
#' @param x,y Numeric vectors.
#' @param r,s Orders.
#' @param central If TRUE, moments about the means.
#' @export
moment_rs <- function(x, y, r, s, central = FALSE) {
  if (central) { x <- x - mean(x); y <- y - mean(y) }
  mean(x^r * y^s)
}

#' Population covariance
#' @param x,y Numeric vectors.
#' @export
cov_pop <- function(x, y) mean(x * y) - mean(x) * mean(y)

#' Linear correlation coefficient (population formulas)
#' @param x,y Numeric vectors.
#' @export
cor_pop <- function(x, y) cov_pop(x, y) / (sd_pop(x) * sd_pop(y))

#' Simple linear regression by hand: y = a + b x
#' @param x,y Numeric vectors.
#' @return List with a, b, r, r2, sse, mse, fitted values and residuals.
#' @examples
#' simple_regression(c(1, 2, 3, 4, 5), c(1, 1, 2, 2, 4))
#' @export
simple_regression <- function(x, y) {
  b <- cov_pop(x, y) / variance_pop(x)
  a <- mean(y) - b * mean(x)
  fit <- a + b * x
  res <- y - fit
  sse <- sum(res^2)
  list(a = a, b = b, r = cor_pop(x, y), r2 = cor_pop(x, y)^2,
       sse = sse, mse = sse / length(y), fitted = fit, residuals = res,
       equation = sprintf("y = %.4f + %.4f x", a, b))
}

#' Fit a model after a variable transformation
#' @description Types: linear (a + b x), exponential (a e^(b x)), power (a x^b),
#' hyperbolic (a + b / x), logarithmic (a + b ln x), parabola (a + b x + c x^2).
#' SSE and R^2 are computed on the ORIGINAL y scale so models can be compared.
#' @param x,y Numeric vectors (positive where logs are needed).
#' @param type Model type.
#' @return List with coefficients, fitted, residuals, sse, mse, r2, equation, predict function.
#' @examples
#' fit_model(1:5, c(2.1, 3.9, 8.2, 15.8, 32.5), "exponential")
#' @export
fit_model <- function(x, y, type = c("linear", "exponential", "power",
                                      "hyperbolic", "logarithmic", "parabola")) {
  type <- match.arg(type)
  d <- data.frame(x = x, y = y)
  co <- switch(type,
    linear = { m <- lm(y ~ x, d); c(a = unname(coef(m)[1]), b = unname(coef(m)[2])) },
    exponential = { m <- lm(log(y) ~ x, d); c(a = exp(unname(coef(m)[1])), b = unname(coef(m)[2])) },
    power = { m <- lm(log(y) ~ log(x), d); c(a = exp(unname(coef(m)[1])), b = unname(coef(m)[2])) },
    hyperbolic = { d$u <- 1 / d$x; m <- lm(y ~ u, d); c(a = unname(coef(m)[1]), b = unname(coef(m)[2])) },
    logarithmic = { d$u <- log(d$x); m <- lm(y ~ u, d); c(a = unname(coef(m)[1]), b = unname(coef(m)[2])) },
    parabola = { d$x2 <- d$x^2; m <- lm(y ~ x + x2, d)
                 c(a = unname(coef(m)[1]), b = unname(coef(m)[2]), c = unname(coef(m)[3])) })
  pred <- function(z) {
    switch(type,
      linear = co[["a"]] + co[["b"]] * z,
      exponential = co[["a"]] * exp(co[["b"]] * z),
      power = co[["a"]] * z^co[["b"]],
      hyperbolic = co[["a"]] + co[["b"]] / z,
      logarithmic = co[["a"]] + co[["b"]] * log(z),
      parabola = co[["a"]] + co[["b"]] * z + co[["c"]] * z^2)
  }
  fit <- pred(x)
  res <- y - fit
  sse <- sum(res^2)
  sst <- sum((y - mean(y))^2)
  list(type = type, coefficients = co, fitted = fit, residuals = res,
       sse = sse, mse = sse / length(y), r2 = 1 - sse / sst,
       equation = paste0(type, ": ", paste(names(co), round(co, 4), sep = "=", collapse = ", ")),
       predict = pred)
}

#' Compare all transformation models (SSE, MSE, R^2)
#' @param x,y Numeric vectors.
#' @param types Models to compare.
#' @return Data frame ordered by decreasing R^2.
#' @export
compare_fits <- function(x, y, types = c("linear", "exponential", "power",
                                          "hyperbolic", "logarithmic", "parabola")) {
  rows <- lapply(types, function(tp) {
    f <- tryCatch(fit_model(x, y, tp), error = function(e) NULL)
    if (is.null(f) || any(!is.finite(f$fitted))) return(NULL)
    data.frame(model = tp, sse = f$sse, mse = f$mse, r2 = f$r2)
  })
  out <- do.call(rbind, rows)
  out[order(-out$r2), ]
}

#' Multiple linear regression summary
#' @param formula A formula, e.g. y ~ x1 + x2.
#' @param data A data frame.
#' @return List with the lm object, coefficients, R^2, adjusted R^2 and residual sd.
#' @export
multiple_regression <- function(formula, data) {
  m <- lm(formula, data = data)
  s <- summary(m)
  list(model = m, coefficients = coef(m), r2 = s$r.squared,
       adj_r2 = s$adj.r.squared, sigma = s$sigma)
}
