#' Computes the exact expectation of the number of true and false rejections for the BH procedure.
#'
#' @param K0 Number of null baskets.
#' @param K1 Number of alternative baskets.
#' @param alpha Threshold parameter for BH.
#' @param F1 The distribution of p-values under the alternative.
#' @param F0 The distribution of p-values under the null..
#'
#' @return A list of the expected number of true and false rejections.
#' @export
exact_bh_expectations <- function(K1, K0, alpha, F1, F0) {
  K <- K1 + K0

  # BH critical values: c_vals[1] = 0, c_vals[2] = alpha/K, ..., c_vals[K+1] = alpha
  c_vals <- (0:K) * alpha / K


  # 1. Initialize joint distribution at maximum threshold c_K using F0 and F1
  prob0 <- dbinom(0:K0, size = K0, prob = F0(c_vals[K + 1]))
  prob1 <- dbinom(0:K1, size = K1, prob = F1(c_vals[K + 1]))

  P <- outer(prob0, prob1)

  n0_grid <- row(P) - 1
  n1_grid <- col(P) - 1

  expected_S <- 0.0
  expected_V <- 0.0

  # 2. Process backwards from K down to 1
  for (k in K:1) {

    # A. Absorption
    absorb_mask <- (P > 0) & ((n0_grid + n1_grid) >= k)
    if (any(absorb_mask)) {
      expected_S <- expected_S + sum(n1_grid[absorb_mask] * P[absorb_mask])
      expected_V <- expected_V + sum(n0_grid[absorb_mask] * P[absorb_mask])
      P[absorb_mask] <- 0.0
    }

    if (k == 1) {
      break
    }

    # B. Backward Transition: Thinning based on F0 and F1
    prob_upper0 <- F0(c_vals[k + 1])
    p0 <- if (prob_upper0 > 0) F0(c_vals[k]) / prob_upper0 else 0.0

    prob_upper1 <- F1(c_vals[k + 1])
    p1 <- if (prob_upper1 > 0) F1(c_vals[k]) / prob_upper1 else 0.0

    # Thin the nulls (n0 -> m0)
    Temp <- matrix(0, nrow = K0 + 1, ncol = K1 + 1)
    for (n0 in 0:K0) {
      if (any(P[n0 + 1, ] > 0)) {
        trans0 <- dbinom(0:n0, size = n0, prob = p0)
        for (m0 in 0:n0) {
          Temp[m0 + 1, ] <- Temp[m0 + 1, ] + P[n0 + 1, ] * trans0[m0 + 1]
        }
      }
    }

    # Thin the alternatives (n1 -> m1)
    NextP <- matrix(0, nrow = K0 + 1, ncol = K1 + 1)
    for (n1 in 0:K1) {
      if (any(Temp[, n1 + 1] > 0)) {
        trans1 <- dbinom(0:n1, size = n1, prob = p1)
        for (m1 in 0:n1) {
          NextP[, m1 + 1] <- NextP[, m1 + 1] + Temp[, n1 + 1] * trans1[m1 + 1]
        }
      }
    }

    P <- NextP
  }

  return(list(
    expected_true_rejections = expected_S,
    expected_false_rejections = expected_V
  ))
}

#' Computes the probability of false rejections for the BH procedure.
#'
#' @param n1 Stage 1 sample size.
#' @param r1 Stage 1 futility threshold.
#' @param n2 Stage 2 sample size.
#' @param p0 Response probability under the null.
#' @param pa Response probability under the alternative.
#' @param K Total number of null baskets.
#' @param K1 Number of alternative baskets.
#' @param alpha Threshold parameter for BH.
#'
#' @return The probability of false rejections.

#' @export
bh_tpr=function(n1,r1,n2,p0,pa,K,alpha,K1)
{

  pmf0=pmf_two_stage(n1, r1, n2, p0)
  pmf1=pmf_two_stage(n1, r1, n2, pa)

  p_val <- Vectorize(function(x) {
    if (x <= r1) prob_reject=1-pbinom(x-1,n1,p0) else {
      prob_reject <- 0
      for (x1 in (r1 + 1):n1) {
        p_x1 <- dbinom(x1, n1, p0)
        p_x2_needed <- 1 - pbinom(x - x1-1, n2, p0)
        prob_reject <- prob_reject + (p_x1 * p_x2_needed)
      } }
    return(prob_reject)
  })

  p_vals=p_val(0:(n1+n2))

  F_0=Vectorize(function(t)
  {
    ind_temp=which(p_vals<=t)
    sum(pmf0[ind_temp])
  })


  F_1=Vectorize(function(t)
  {
    ind_temp=which(p_vals<=t)
    sum(pmf1[ind_temp])
  })

  #plot(F_0,0,1)
  #abline(0,1)



  out=exact_bh_expectations(K1=K1, K0=K-K1, alpha=alpha, F1=F_1, F0=F_0)
  out$expected_true_rejections/K1
}

