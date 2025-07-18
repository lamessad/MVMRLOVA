#' Example Dataset for MVMRLOVA
#'
#' This dataset contains simulated data for multivariable Mendelian Randomization (MVMR) analysis.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{betaY}{A numeric vector of length \eqn{n}, where \eqn{n} is the number of instrumental variables (IVs), containing the effect sizes of IVs on the outcome}
#'   \item{betaYse}{A numeric vector of length \eqn{n}, containing the standard errors of the outcome effect sizes}
#'   \item{betaX}{A numeric matrix of dimension \eqn{n \times p}, where \eqn{p} is the number of exposures; each row corresponds to an IV and each column to an exposure's effect size}
#'   \item{betaXse}{A numeric matrix of dimension \eqn{n \times p}, containing the standard errors corresponding to \code{betaX}}
#'   \item{ny}{A single numeric value giving the sample size of the outcome GWAS}
#' }
#' @source The data was simulated/generated for illustrative purposes.
"mvmrlova_dat"

