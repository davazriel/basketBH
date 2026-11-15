# ------------------------------------------------------------
# pmf_two_stage(n1, r1, n2, p)
# Returns: named vector of probabilities for X (two-stage experiment)
# ------------------------------------------------------------
pmf_two_stage <- function(n1, r1, n2, p) {
  xs <- 0:(n1 + n2)
  pmf <- rep(0, length(xs))

  # Stage 1 distribution
  for (x1 in 0:n1) {
    px1 <- dbinom(x1, size = n1, prob = p)

    if (x1 <= r1) {
      # Stop immediately
      pmf[x1 + 1] <- pmf[x1 + 1] + px1
    } else {
      # Continue to stage 2, X = x1 + X2
      for (x2 in 0:n2) {
        pmf[x1 + x2 + 1] <- pmf[x1 + x2 + 1] +
          px1 * dbinom(x2, size = n2, prob = p)
      }
    }
  }
  names(pmf) <- xs
  pmf
}


expected_sample_size <- function(K, r1, n1, n2, p0) {
  # Probability of an inactive arm passing stage 1
  p_pass <- 1 - pbinom(r1, n1, p0)
  # Expected total sample size
  en0 <- K * n1 + K * n2 * p_pass
  return(en0)
}

