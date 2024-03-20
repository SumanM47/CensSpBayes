#'simdat: A simulated spatial data with left-censoring
#'
#' Simulated Dataset to use the CensSpBayes function on
#'
#'
#' @format ## `simdat`
#' A data frame with 10,000 rows and 4 columns:
#' \describe{
#'   \item{Y}{The response, censored whenever this is smaller than cutoff.Y}
#'   \item{cutoff.Y}{The limit below which any data is censored}
#'   \item{x, y}{x and y locations of the observations}
#' }
"simdat"
