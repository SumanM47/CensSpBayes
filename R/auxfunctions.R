#' @name transform
#' @title Internal list to compute certain transformations
#'
#' @param x a real vector, within upper and lower limits for logit and any real value for inv.logit
#' @param lower lower limit for the logit and expit transform, defaults to 0
#' @param upper upper limit for the logit and expit transform, defaults to 1
#'
#' @return list of function
#' @noRd

transform <- list(
  logit = function(x, lower = 0, upper = 1) {
    x <- (x - lower) / (upper - lower)
    return(log(x / (1 - x)))
  },
  inv.logit = function(x, lower = 0, upper = 1) {
    p <- exp(x) / (1 + exp(x))
    p <- p * (upper - lower) + lower
    return(p)
  }
)

#' @name mhupdate
#' @title Internal function to correct proposal variances for the MH updates to keep the acceptance ratio between 0.3 and 0.5
#'
#' @param acc integer, number of times the proposed variable has been accepted
#' @param att integer, number of attempts made
#' @param mh positive scalars, proposal variances
#' @param nattempts integer, minimum number of attempts that need to be made before updating the proposal variances
#' @param lower number between 0 and 1, the ratio by which to decrease the proposal variance should it be necessary
#' @param upper number between 1 and infinity, the ratio by which to increase the proposal variance should it be necessary
#'
#' @return a list of the following variables: acc (integer vector, set to 0 if enough number of attempts were made), att (integer vector, set to 0 if enough number of attempts were made), mh (positive scalars, possibly updated proposal variances for the MH update)
#'
#' @noRd

mhupdate <- function(acc, att, mh, nattempts = 50, lower = 0.8, higher = 1.2) {
  acc.rate     <- acc / att
  these.update <- att > nattempts
  these.low    <- (acc.rate < 0.30) & these.update
  these.high   <- (acc.rate > 0.50) & these.update

  mh[these.low]  <- mh[these.low] * lower
  mh[these.high] <- mh[these.high] * higher

  acc[these.update] <- 0
  att[these.update] <- 0

  results <- list(acc = acc, att = att, mh = mh)
  return(results)
}
