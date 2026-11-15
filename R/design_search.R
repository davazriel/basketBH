#' Search for BH Basket Designs
#'
#' Searches for two-stage basket trial designs that satisfy
#' a target power requirement under the BH multiple-testing
#' procedure.
#'
#' @param K Number of baskets.
#' @param p0 Response probability under the null.
#' @param pa Response probability under the alternative.
#' @param target_alpha Target FDR level.
#' @param target_power Desired power.
#' @param min_n Minimum total sample size.
#' @param max_n Maximum total sample size.
#' @param progress Logical; print progress updates.
#' @param progress_step Percentage increment for reporting.
#'
#' @return A list containing the optimal and minimax designs.
#'
#' @export
# Main function
design_bh <- function(K, p0, pa,
                     target_alpha = 0.1,
                     target_power = 0.80,
                     min_n = 10,
                     max_n = 40,
                     progress = TRUE,
                     progress_step = 1) {

  # Data frame to store valid designs
  valid_designs <- data.frame()

  # ----------------------------------------------------------
  # Set up progress reporting
  #
  # For fixed n and n1, r1 ranges from 0 to n1, giving
  # n1 + 1 computations.
  # ----------------------------------------------------------
  total_computations <- sum(
    vapply(min_n:max_n, function(n) {
      sum((1:(n - 1)) + 1)
    }, numeric(1))
  )

  completed_computations <- 0

  # Percentage at which the next update will be printed
  next_report <- progress_step

  if (progress) {
    if (progress_step <= 0 || progress_step > 100) {
      stop("progress_step must be between 0 and 100.")
    }

    cat(sprintf(
      "Starting search over %d candidate designs...\n",
      total_computations
    ))
  }


  # Loop over total sample size n
  for (n in min_n:max_n) {
    # Loop over stage 1 sample size
    for (n1 in 1:(n - 1)) {

      n2 <- n - n1

      # Loop over stage 1 threshold
      for (r1 in 0:n1) {

        # Compute power
        expected_power <- bh_tpr(
          n1, r1, n2, p0, pa, K, target_alpha, 1)

        # Update the progress counter
        completed_computations <- completed_computations + 1

        if (progress) {
          percentage_completed <-
            100 * completed_computations / total_computations

          if (percentage_completed >= next_report) {
            cat(sprintf(
              "\rProgress: %6.2f%% (%d of %d computations)",
              percentage_completed,
              completed_computations,
              total_computations
            ))

            flush.console()

            # Move next_report beyond the current percentage.
            # This also handles progress_step values below 1.
            while (next_report <= percentage_completed) {
              next_report <- next_report + progress_step
            }
          }
        }

        # Check whether the power constraint is met
        if (expected_power >= target_power) {

          # Calculate expected sample size per basket under p0
          en0 <- expected_sample_size(K, r1, n1, n2, p0) / K

          # Save the valid design
          valid_designs <- rbind(
            valid_designs,
            data.frame(
              n = n,
              n1 = n1,
              n2 = n2,
              r1 = r1,
              EN0 = en0
            )
          )
        }
      }
    }
  }

  if (progress) {
    cat(sprintf(
      "\rProgress: 100.00%% (%d of %d computations)\n",
      total_computations,
      total_computations
    ))
  }

  # Extract the winners
  if (nrow(valid_designs) > 0) {

    # Optimal design: minimum EN0
    ind.optimal <- which.min(valid_designs$EN0)
    optimal <- valid_designs[ind.optimal, , drop = FALSE]
    rownames(optimal) <- NULL

    # Minimax design: minimum n, and then minimum EN0
    minimax_subset <-
      valid_designs[valid_designs$n == min(valid_designs$n), ]

    ind.minimax <- which.min(minimax_subset$EN0)
    minimax <- minimax_subset[ind.minimax, , drop = FALSE]
    rownames(minimax) <- NULL

    cat("\n--- OPTIMAL DESIGN ---\n")
    print(optimal, row.names = FALSE)


    cat("\n--- MINIMAX DESIGN ---\n")
    print(minimax, row.names = FALSE)

    return(list(
      optimal = optimal,
      minimax = minimax
    ))

  } else {
    cat("No designs found within the limits.\n")
    return(NULL)
  }
}
