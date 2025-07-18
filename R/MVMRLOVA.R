
#' MultiVariable Mendelian Randomization based on Latent Outcome Variable Approach(MVMRLOVA)
#'
#' Identifies vertical pleiotropy variants and estimates multivariable causal effects using GWAS summary statistics.
#'
#' @param betaY A numeric vector of standardized effect sizes (beta) from the outcome GWAS.
#' @param betaX A numeric matrix of standardized effect sizes (beta) from the exposures GWAS. Each column corresponds to a different exposure.
#' @param betaYse A numeric vector of standard errors for `betaY`.
#' @param betaXse A numeric matrix of standard errors for `betaX`.
#' @param ny sample size of outcome GWAS
#' @param gwas_p p-value threshold for the association between variants and the latent phenotype (the direct genetic effect). Default: 5e-2.
#' @param gwas_p2 p-value threshold for exposure GWAS used in instrument variable selection. Default: 5e-8
#' @param permutn number of permutations, Default: 0.
#' @param log_file file name for saving the iterations to check convergence. Default: NULL.
#'
#' @details
#' This function iteratively updates the set of vertical pleiotropy instruments (Valid) by identifying SNPs that pass significance criteria based on latent outcome variable GWAS statistics.
#' It assumes that both SNP effect sizes and the outcome trait are standardized.
#'
#' The function uses the `MendelianRandomization` package to perform multivariable inverse-variance weighted (IVW) MR.
#'
#' @return A list with the following elements:
#' \item{CausEst}{Estimate of causal effect.}
#' \item{CausEstSE}{Standard error of causal effect estimate.}
#' \item{CausEstP}{P-value from the z test for the causal effect.}
#' \item{IVs}{p-values for the direct causal effect of all Instrumental Variants (IV) on the outcome.}
#' \item{Valid}{Index of valid  instrumental variants.}
#' \item{sig_v}{the 5th percentile of the permuted p-values .}
#' \item{corrected_p-value}{corrected p-value based on the permutation distribution.}
#'@return An object of class "mvmrlova" with summary output.
#'
#' @keywords Mendelian randomization.
#' @import MendelianRandomization
#' @importFrom stats cor ecdf lm pchisq pt quantile
#'
#' @export
#'
#'
MVMRLOVA = function(betaY, betaX, betaYse, betaXse, ny, gwas_p , gwas_p2, permutn = 0, log_file = NULL) {


  log_message <- function(message, log_file) {
    if (!is.null(log_file)) {
      cat(message, file = log_file, append = TRUE, sep = "\n")
    }
  }

  MVMR_LOVA = function(betaY, betaX, betaYse, betaXse, ny,  gwas_p, gwas_p2, log_file = NULL) {
    sv1 = seq(1, length(betaY))


    tmp = cbind(betaY, betaX)
    pcor = cor(tmp)
    nexp = ncol(betaX)

    pX = pchisq((betaX / betaXse)^2, 1, lower.tail = FALSE)

    yi=0;tau_t=array(-9999,nexp);tau_t2=array(0,nexp)


    while (!all(tau_t == tau_t2) && yi <= 10) {

      yi = yi + 1
      tau_t = tau_t2

      gwas5 = matrix(0, length(betaY), 4)
      gwas5[,1] = betaY - betaX %*% tau_t

      cov_tmp = 0
      for (i in 1:nexp) {
        cov_tmp = cov_tmp - 2 * tau_t[i] * pcor[1, 1 + i]
        for (j in 1:nexp) {
          if (i != j) {
            cov_tmp = cov_tmp + 2 * tau_t[i] * tau_t[j] * pcor[1 + i, 1 + j]
          }
        }
      }

      gwas5[,2] = sqrt(betaYse^2 + sum(tau_t^2 / ny) - cov_tmp / ny)
      gwas5[,3] = gwas5[,1] / gwas5[,2]
      gwas5[,4] = pchisq(gwas5[,3]^2, 1, lower.tail = FALSE)

      svv1 = sv1[gwas5[,4] > gwas_p]
      svv2 = which(apply(pX, 1, function(row) any(row < gwas_p2)))
      svv = intersect(svv1, svv2)

      mv_ivw <- MendelianRandomization::mr_mvinput(bx = betaX[svv,], bxse = betaXse[svv,],
                           by = betaY[svv], byse = betaYse[svv])
      ivw <- MendelianRandomization::mr_mvivw(mv_ivw)
      tau_t2 = ivw$Estimate


    }

    mv_ivw <- MendelianRandomization::mr_mvinput(bx = betaX[svv,], bxse = betaXse[svv,],
                         by = betaY[svv], byse = betaYse[svv])
    ivw <- MendelianRandomization::mr_mvivw(mv_ivw)

    return(list(
      CausEst = ivw$Estimate,
      CausEstSE = ivw$StdError,
      CausEstP = ivw$Pvalue,
      IVs = gwas5[,4],
      Valid = svv,
      pcor=pcor
    ))
  }

  mvmr_lova_result <- MVMR_LOVA(betaY, betaX, betaYse, betaXse,  ny, gwas_p, gwas_p2, log_file)

  exposure_names <- colnames(betaX)

  if (permutn > 0) {
    permutp = matrix(NA, permutn, ncol(betaX))

    for (i in 1:permutn) {
      sv2 = sample(seq_len(length(betaY)))
      perm_result = suppressWarnings(MVMR_LOVA(betaY[sv2], betaX, betaYse[sv2], betaXse,  ny,  gwas_p, gwas_p2))
      permutp[i, ] = perm_result$CausEstP
    }

    permutt = apply(permutp, 2, function(pnull) quantile(pnull, probs = 0.05, na.rm = TRUE))
    corrected_p = mapply(function(pval, null) ecdf(null)(pval),
                         pval = mvmr_lova_result$CausEstP,
                         null = as.data.frame(permutp))

    if (any(1 / mvmr_lova_result$CausEstP > permutn)) {
      warning("To get a more precise p-value, it is recommended to increase the number of permutations.")
      log_message("Recommendation: Increase permutation number.", log_file)
    }

    result <- list(
      CausEst = mvmr_lova_result$CausEst,
      CausEstSE = mvmr_lova_result$CausEstSE,
      CausEstP = mvmr_lova_result$CausEstP,
      IVs = mvmr_lova_result$IVs,
      Valid = mvmr_lova_result$Valid,
      sig_v = permutt,
      corrected_p = corrected_p,
      permutn = permutn,
      exposure_names = exposure_names
    )
  } else {
    result <- list(
      CausEst = mvmr_lova_result$CausEst,
      CausEstSE = mvmr_lova_result$CausEstSE,
      CausEstP = mvmr_lova_result$CausEstP,
      IVs = mvmr_lova_result$IVs,
      Valid = mvmr_lova_result$Valid,
      sig_v = NULL,
      corrected_p = NULL,
      permutn = 0,
      exposure_names = exposure_names
    )
  }

  class(result) <- "mvmrlova"
  return(result)
}
