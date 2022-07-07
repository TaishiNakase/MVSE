#' MVSE: A package for estimating (M)osquito-borne (V)iral (S)uitability (E)stimator
#'
#' The MVSE package provides computational methods to estimate an index of 'transmission potential' for mosquito-borne viruses. MVSE allows for
#' parameterization of particular viruses, host and mosquito-species. The methods
#' offered are climate-driven and can therefore be location specific.
#'
#' For more information, please refer to the manual.
#' 
#' @docType package
#' @author Taishi Nakase
#' @import Rcpp
#' @importFrom  Rcpp evalCpp
#' @useDynLib MVSE
#' @exportPattern "^[[:alpha:]]+"
#' @include mvsemodel.R mvsefit.R utils.R
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
