#' Print method for MVMRLOVA result
#'
#' @param x An object of class \code{mvmrlova}
#' @param ... Additional arguments (ignored)
#'
#' @export
#'
print.mvmrlova <- function(x, ...) {
  cat("\nMultivariable MR via Latent Outcome Variable Approach (MVMRLOVA)\n")
  cat("--------------------------------------------------------------------------\n")
  cat("Number of SNPs considered     :", length(x$IVs), "\n")
  cat("Number of valid instruments   :", length(x$Valid), "\n")

  if (!is.null(x$sig_v)) {
    cat("Number of permutations        :", x$permutn, "\n")
  }

  cat("--------------------------------------------------------------------------\n")

  est <- x$CausEst
  se <- x$CausEstSE
  pval <- x$CausEstP
  lower <- est - 1.96 * se
  upper <- est + 1.96 * se

  exposure_labels <- if (!is.null(x$exposure_names)) x$exposure_names else paste0("Exposure_", seq_along(est))

  result_table <- data.frame(
    Exposure = exposure_labels,
    Estimate = round(est, 3),
    `Std. Error` = round(se, 3),
    `95% CI` = paste0("(", round(lower, 3), ", ", round(upper, 3), ")"),
    `p-value` = formatC(pval, format = "e", digits = 2),
    stringsAsFactors = FALSE
  )

  if (!is.null(x$corrected_p)) {
    result_table$`Corrected p` <- formatC(x$corrected_p, format = "e", digits = 2)
  }

  print(result_table, row.names = FALSE)
}
