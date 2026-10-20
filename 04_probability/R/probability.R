#' Classical (Laplace) probability
#' @param favourable Number of favourable cases (can be choose(...) * choose(...)).
#' @param possible Number of possible cases.
#' @examples
#' classical_prob(choose(3, 1) * choose(17, 3), choose(20, 4))  # 0.42
#' @export
classical_prob <- function(favourable, possible) favourable / possible

#' Probability of the complement
#' @param pA P(A).
#' @export
prob_complement <- function(pA) 1 - pA

#' Probability of a union: P(A) + P(B) - P(A and B)
#' @param pA P(A).
#' @param pB P(B).
#' @param pAB P(A and B); 0 if mutually exclusive.
#' @export
prob_union <- function(pA, pB, pAB = 0) pA + pB - pAB

#' Conditional probability P(A|B) = P(A and B) / P(B)
#' @param pAB P(A and B).
#' @param pB P(B), must be > 0.
#' @export
cond_prob <- function(pAB, pB) {
  if (pB <= 0) stop("P(B) must be positive.")
  pAB / pB
}

#' Check independence of two events
#' @param pA P(A).
#' @param pB P(B).
#' @param pAB P(A and B).
#' @param tol Numerical tolerance.
#' @return TRUE if P(A and B) = P(A) P(B).
#' @export
is_independent <- function(pA, pB, pAB, tol = 1e-9) abs(pAB - pA * pB) < tol

#' Theorem of total probability
#' @param prior Vector P(A_i) of a complete system (sums to 1).
#' @param likelihood Vector P(B|A_i).
#' @examples
#' total_probability(c(0.2, 0.8), c(0.9, 0.4))  # 0.5
#' @export
total_probability <- function(prior, likelihood) {
  if (abs(sum(prior) - 1) > 1e-9) warning("Priors do not sum to 1.")
  sum(prior * likelihood)
}

#' Bayes' theorem: posterior probabilities P(A_i|B)
#' @param prior Vector P(A_i).
#' @param likelihood Vector P(B|A_i).
#' @return Vector of posteriors (sums to 1).
#' @examples
#' bayes(c(0.2, 0.8), c(0.9, 0.4))  # 0.36 0.64
#' @export
bayes <- function(prior, likelihood) {
  joint <- prior * likelihood
  joint / sum(joint)
}

#' Diagnostic test measures
#' @param prevalence P(disease).
#' @param sensitivity P(+|disease).
#' @param specificity P(-|no disease).
#' @return List with P(+), PPV = P(disease|+) and NPV = P(no disease|-).
#' @examples
#' diagnostic_test(0.10, 0.95, 0.90)  # PPV = 0.5135
#' @export
diagnostic_test <- function(prevalence, sensitivity, specificity) {
  p_pos <- prevalence * sensitivity + (1 - prevalence) * (1 - specificity)
  p_neg <- 1 - p_pos
  list(p_positive = p_pos,
       ppv = prevalence * sensitivity / p_pos,
       npv = (1 - prevalence) * specificity / p_neg)
}

#' Number of independent repetitions needed to reach a target probability
#' @param p Probability of the event in one trial.
#' @param target Wanted probability of at least one success (default 0.5).
#' @examples
#' repetitions_needed(1/36)  # 25 rolls of two dice for at least one double six
#' @export
repetitions_needed <- function(p, target = 0.5) {
  ceiling(log(1 - target) / log(1 - p))
}
