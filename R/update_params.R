#' @name cormat.inv.update.inla
#' @title Internal function to compute the inverse and log-determinant of the approximate correlation matrix using SPDE
#'
#' @param rho positive scalar, spatial range parameter
#' @param c.mat the mass matrix for the SPDE approximation
#' @param g1.mat the stiffness matrix for the SPDE approximation
#' @param g2.mat G2 matrix created by the mass and stiffness matrices
#' @param alpha order of SPDE, defaults to 2, the only supported value at the time
#'
#' @import spam
#'
#' @return list of quantities: cormat.inv(approximated inverse of the correlation matrix), cormat.logdet(approximated log-determinant for the correlation matrix)
#' @noRd

cormat.inv.update.inla <- function(rho, c.mat, g1.mat, g2.mat, alpha = 2){
  cormat.inv <- (1 / rho)^4 * c.mat + 2 * (1 / rho)^2 * g1.mat + g2.mat
  tau <- rho^2 / (4 * pi)
  cormat.inv <- tau * cormat.inv
  cormat.logdet <- -2 * sum(log(diag(spam::chol(cormat.inv))))
  list(cormat.inv = cormat.inv, cormat.logdet = cormat.logdet)}

#' @name theta.latent.update
#' @title Internal function to update the coefficient vectors and the latent parameters simultaneously
#'
#' @param nq integer, number of columns of the design matrix
#' @param nmesh integer, the size of the mesh used for the SPDE solver
#' @param crossprod.X X'X matrix, X being the design matrix
#' @param crossprod.A A'A matrix, A being the index matrix for SPDE
#' @param crossprod.A.X A'X matrix
#' @param crossprod.X.A X'A matrix
#' @param crossprod.Xy X'y vector, y being the spatial observations
#' @param crossprod.Ay A'y vector
#' @param tau positive scalar, total precision
#' @param cormat.inv inverted correlation matrix
#' @param r positive scalar, the ratio of partial sill to sill
#' @param sd.theta prior standard deviation for theta
#'
#' @import spam
#' @import stats
#'
#' @return updated vector of size nq+nmesh with the first nq many being coefficients, the latter nmesh being the latent parameters
#' @noRd

theta.latent.update <- function(nq, nmesh, crossprod.X, crossprod.A,
                                crossprod.A.X, crossprod.X.A,
                                crossprod.Xy, crossprod.Ay,
                                tau, cormat.inv, r, sd.theta){

  latent.cov.inv.11 <- crossprod.X / (1 - r) + diag(nq) / sd.theta^2
  latent.cov.inv.12 <- crossprod.X.A / (1 - r)
  latent.cov.inv.21 <- crossprod.A.X / (1 - r)
  latent.cov.inv.22 <- crossprod.A / (1 - r) + cormat.inv / r

  latent.mean.part1 <- crossprod.Xy / (1 - r)
  latent.mean.part2 <- crossprod.Ay / (1 - r)

  latent.cov.inv <- rbind(cbind(latent.cov.inv.11, latent.cov.inv.12),
                          cbind(latent.cov.inv.21, latent.cov.inv.22))
  latent.mean.part <- c(latent.mean.part1, latent.mean.part2)

  chol.latent.cov.inv <- spam::chol(latent.cov.inv)
  tchol.latent.cov.inv <- t(chol.latent.cov.inv)
  omega <- spam::forwardsolve(tchol.latent.cov.inv, latent.mean.part)
  mm <- spam::backsolve(chol.latent.cov.inv, omega)
  zz <- rnorm(nq + nmesh)
  vv <- spam::backsolve(chol.latent.cov.inv, zz)
  theta.latent <- mm + vv / sqrt(tau)
  theta.latent}

#' @name tau.update
#' @title Internal function to update the precision parameter
#'
#' @param ns integer, number of observed locations
#' @param nq integer, number of columns for the design matrix
#' @param nmesh integer, size of the mesh for SPDE solver
#' @param cur.rss positive scalar, current residual sum of squares
#' @param cur.ss.theta positive scalar, current sum of squares due to the coefficients
#' @param cur.ss.latent positive scalar, current sum of squares for the latent parameters
#' @param r proportion of partial sill to sill
#' @param tau.a hyperparameter for the inverse gamma prior for tau
#' @param tau.b hyperparameter for the inverse gamma prior for tau
#'
#' @import stats
#'
#' @return positive scalar, new precision parameter estimate
#' @noRd

tau.update <- function(ns, nq, nmesh, cur.rss, cur.ss.theta, cur.ss.latent, r, tau.a, tau.b){
  tau <- rgamma(1, shape = tau.a + (ns + nq + nmesh) / 2,
                rate = tau.b + (cur.rss / (1 - r) + cur.ss.theta + cur.ss.latent / r) / 2)
  tau}


#' @name rho.update
#' @title Internal function to update the spatial range parameter
#'
#' @param rho positive scalar, current value of the spatial range parameter
#' @param c.mat the mass matrix for the SPDE approximation
#' @param g1.mat the stiffness matrix for the SPDE approximation
#' @param g2.mat G2 matrix created by the mass and stiffness matrices
#' @param latent vector of latent parameters
#' @param tau positive scalar, the precision
#' @param r proportion of partial sill to sill
#' @param cormat.inv inverse-correlation matrix
#' @param cormat.logdet log-determinant of the correlation matrix
#' @param cur.ss.latent positive scalar, current sum of squares due to latent parameters
#' @param rho.upper hyperparameter corresponding to the upper limit for the uniform prior on rho, defaults to infinity
#' @param att.rho integer, number of attempts to update rho
#' @param acc.rho integer, number of times proposed rho value has been updated
#' @param mh.rho positive scalar, proposal variance for rho update
#'
#' @import stats
#' @import spam
#'
#' @return list of quantities: rho(updated spatial range parameter), cormat.inv(updated inverse correlation matrix), cormat.logdet(updated log-determinant for the correlation matrix), cur.ss.latent (updated current sum of squares due to latent parameters), att.rho (updated number of attempts made to update rho), acc.rho (updated number of times proposed rho was accepted)
#' @noRd



rho.update <- function(rho, c.mat, g1.mat, g2.mat, latent,
                       tau, r, cormat.inv, cormat.logdet, cur.ss.latent,
                       rho.upper = Inf, att.rho, acc.rho, mh.rho){

  att.rho <- att.rho + 1

  rho.star <- transform$logit(rho, lower = 0, upper = rho.upper)
  can.rho.star <- rnorm(1, rho.star, mh.rho)
  can.rho <- transform$inv.logit(can.rho.star, lower = 0, upper = rho.upper)

  can.cormat.details <- cormat.inv.update.inla(can.rho, c.mat, g1.mat, g2.mat)
  can.cormat.inv <- can.cormat.details$cormat.inv
  can.cormat.logdet <- can.cormat.details$cormat.logdet

  can.ss.latent <- sum((can.cormat.inv %*% latent) * latent)

  ratio <- -0.5 * tau * (can.ss.latent - cur.ss.latent) / r - 0.5 *
    (can.cormat.logdet - cormat.logdet) +
    log(can.rho - 0) + log(rho.upper - can.rho) -
    log(rho - 0) - log(rho.upper - rho)

  if(log(runif(1)) < ratio){
    rho <- can.rho
    cormat.inv <- can.cormat.inv
    cormat.logdet <- can.cormat.logdet
    acc.rho <- acc.rho + 1
    cur.ss.latent <- can.ss.latent}

  results <- list(rho = rho, cormat.inv = cormat.inv,
                  cormat.logdet = cormat.logdet,
                  cur.ss.latent = cur.ss.latent,
                  att.rho = att.rho,  acc.rho = acc.rho)
  results}

#' @name r.update
#' @title Internal function to update the r parameter
#'
#' @param ns integer, number of spatial locations
#' @param nmesh integer, size of the mesh
#' @param r current value of the r parameter
#' @param tau current value of the tau parameter
#' @param cur.rss current residual sum of squares
#' @param cur.ss.latent current sum of squares due to latent parameters
#' @param att.r integer, number of attempts made to update r parameter
#' @param acc.r integer, number of times r parameter update has been accepted
#' @param mh.r positive scalar, proposal variance for the MH update of the r parameter
#'
#' @import stats
#'
#' @return list of quantities: r(updated value of the r parameter), att.r (updated number of attempts made to update r), acc.r (updated number of times proposed r was accepted)
#' @noRd

r.update <- function(ns, nmesh, r, tau, cur.rss, cur.ss.latent, att.r, acc.r, mh.r){
  att.r <- att.r + 1

  r.star <- transform$logit(r, 1e-4, 0.9999)
  can.r.star <- rnorm(1, r.star, mh.r)
  can.r <- transform$inv.logit(can.r.star, 1e-4, 0.9999)

  ratio1 <- -0.5 * nmesh * (log(can.r) - log(r)) - 0.5 * tau * cur.ss.latent * (1 / can.r - 1 / r)
  ratio2 <- -0.5 * ns * (log(1 - can.r) - log(1 - r)) - 0.5 * tau * cur.rss * (1 / (1 - can.r) - 1 / (1 - r))
  ratio <- ratio1 + ratio2 + log(can.r - 1e-4) + log(0.9999 - can.r) - log(r - 1e-4) - log(0.9999 - r)

  if(log(runif(1)) < ratio){
    r <- can.r
    acc.r <- acc.r + 1}
  results <- list(r = r, att.r = att.r, acc.r = acc.r)
  return(results)}

#' @name fill.censored.y
#' @title Internal function to fill the censored observations with synthetic ones
#'
#' @param y vector of observed variables
#' @param cutoff.y vector of cutoff below which data were censored
#' @param censored.case index vector of whether data were censored or not
#' @param X.theta X%*%theta, the total effect of the covariates
#' @param A.latent A%*%latent, the effect of the latent parameters at the observation locations
#' @param tau precision parameter
#' @param r proportion of partial sill to sill
#'
#' @import stats
#'
#' @return vector of observations with the censored data replaced with synthetic data
#' @noRd

fill.censored.y <- function(y, cutoff.y, censored.cases, X.theta, A.latent, tau, r){
  total.mean <- X.theta + A.latent
  excess <- cutoff.y - total.mean
  total.mean <- total.mean[censored.cases]
  excess <- excess[censored.cases]

  runif.probs <- pnorm(sqrt(tau / (1 - r)) * excess) * runif(length(censored.cases))
  y[censored.cases] <- total.mean + sqrt((1 - r) / tau) * qnorm(runif.probs)
  y}
