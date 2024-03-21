#' @name CensSpBayes
#' @title Function to run MCMC for left censored massive spatial data
#'
#' @description
#' Generate posterior samples for the parameters in the model and predictions of the true process at the prediction locations with pointwise variances for them
#'
#' @usage CensSpBayes(Y, S, X, cutoff.Y, S.pred, X.pred, inla.mats, alpha = 2,
#'                theta.init = NULL, latent.init = NULL, tau.init = NULL,
#'                rho.init = 0.5, r.init = NULL,
#'                jitter=0,
#'                mean.theta = 0, sd.theta = 1e2,
#'                tau.a = 0.1, tau.b = 0.1,
#'                rho.upper = NULL,
#'                iters = 4000, burn = 2000, thin = 5)
#'
#' @param Y vector of observations
#' @param S matrix of observation locations with each row corresponding to one location
#' @param X design matrix for covariates, includes intercept
#' @param cutoff.Y vector of cutoff values below which data are censored
#' @param S.pred matrix of prediction locations, each row corresponds to a location
#' @param X.pred matrix of covariate values at the prediction locations
#' @param inla.mats list of matrices needed for the SPDE approximation. The entries are: c.mat(mass matrix), g1.mat(stiffness matrix), g2.mat(a combination of mass and stiffness matrices), A(index matrix for observation locations), A.pred(index matrix for prediction locations)
#' @param alpha the order of the SPDE, defaults to 2, currently the only accepted value
#' @param theta.init initial values for theta, the coefficients for X
#' @param latent.init initial values for latent parameters
#' @param tau.init initial value for the precision parameter
#' @param rho.init initial value for the spatial range parameter
#' @param r.init initial value for the r parameter, the proportion of partial sill to sill
#' @param jitter numeric variable to indicate whether the initial values are to be used exactly (jitter=0) or to start with a slightly different value (jitter !=0)
#' @param mean.theta hyperparameter corresponding to the mean of the normal prior for theta
#' @param sd.theta hyperparameter corresponding to the sd of the normal prior for theta
#' @param tau.a hyperparameter corresponding to the shape of the inverse gamma prior for tau
#' @param tau.b hyperparameter corresponding to the rate of the inverse gamma prior for tau
#' @param rho.upper hyperparameter corresponding to the upper limit of the uniform prior for rho. Recommend use no more than 0.25*maximum distance of the domain
#' @param iters number of posterior samples to be drawn
#' @param burn number of posterior samples to be discarded as burn-in period samples
#' @param thin thinning interval. The toatl number of iterations is thin*iters
#'
#' @import stats
#' @import spam
#' @importFrom Matrix diag
#'
#' @return a list of posterior samples for theta, tau, rho and r. Additionally, includes the posterior mean and variances of the latent process and the predicted process, and the computation time in minutes
#' @export

CensSpBayes <- function(Y, S, X, cutoff.Y, S.pred, X.pred, inla.mats, alpha = 2,
                    theta.init = NULL, latent.init = NULL, tau.init = NULL,
                    rho.init = 0.5, r.init = NULL,
		                jitter=0,
                    # priors
                    mean.theta = 0, sd.theta = 1e2,
                    tau.a = 0.1, tau.b = 0.1,
                    rho.upper = NULL,
                    # mcmc settings
                    iters = 4000, burn = 2000, thin = 5){

  tick <- proc.time()[3]


  c.mat <- inla.mats$c.mat
  g1.mat <- inla.mats$g1.mat
  g2.mat <- inla.mats$g2.mat
  A <- inla.mats$A
  A.pred <- inla.mats$A.pred

  X <- as.matrix(X)
  X.pred <- as.matrix(X.pred)

  ns <- length(Y)
  nmesh <- ncol(A)
  nq <- ncol(X)
  np <- nrow(S.pred)

  censored.cases <- which(Y <= cutoff.Y)
  Y[censored.cases] <- cutoff.Y[censored.cases]

  if(is.null(theta.init)){
    theta <- c(solve(crossprod(X)+1e-6*diag(nq)) %*% crossprod(X, Y))
  }else{theta <- theta.init}
  X.theta <- c(X %*% theta)

  if(is.null(latent.init)){
    latent <- rep(0, nmesh)
  }else{latent <- latent.init}

  if(is.null(tau.init)){
    tau <- 1 / var(Y - X.theta)
  }else{tau <- tau.init}

  if(is.null(r.init)){
    r <- 0.5
  }else{r <- r.init}

  rho <- rho.init
  if(is.null(rho.upper)){rho.upper <- Inf}


  cormat.details <- cormat.inv.update.inla(rho, c.mat, g1.mat, g2.mat, alpha = alpha)

  cormat.inv <- cormat.details$cormat.inv
  cormat.logdet <- cormat.details$cormat.logdet

  crossprod.X <- crossprod(X)
  crossprod.A <- crossprod(A)
  crossprod.A.X <- crossprod(A, X)
  crossprod.X.A <- crossprod(X, A)
  crossprod.Xy <- as.vector(crossprod(X, Y))
  crossprod.Ay <- as.vector(crossprod(A, Y))


  acc.rho <- att.rho <- mh.rho <- 1
  acc.r <- att.r <- mh.r <- 1

  latent.sum <- rep(0, nmesh)
  latent2.sum <- rep(0, nmesh)

  Y.pred.sum <- rep(0, np)
  Y.pred2.sum <- rep(0, np)

  keepers.theta <- matrix(NA, nrow = iters, ncol = nq)
  keepers.tau <- rep(NA, iters)
  keepers.rho <- rep(NA, iters)
  keepers.r <- rep(NA, iters)

  if(jitter > 0){r <- jitter(r); rho <- jitter(rho); tau <- jitter(tau); theta <- jitter(theta)}

  return.iters <- (burn + 1):iters

  for(iter in 1:iters){for(ttt in 1:thin){
    theta.latent <- theta.latent.update(nq, nmesh, crossprod.X, crossprod.A,
                                        crossprod.A.X, crossprod.X.A,
                                        crossprod.Xy, crossprod.Ay,
                                        tau, cormat.inv, r, sd.theta)
    theta <- theta.latent[1:nq]
    latent <- theta.latent[-(1:nq)]
    X.theta <- as.vector(X %*% theta)
    A.latent <- as.vector(A %*% latent)
    cur.rss <- sum((Y - X.theta - A.latent)^2)
    cur.ss.theta <- sum(theta^2) / sd.theta^2
    cur.ss.latent <- sum((cormat.inv %*% latent) * latent)

    tau <- tau.update(ns, nq, nmesh, cur.rss, cur.ss.theta, cur.ss.latent, r, tau.a, tau.b)

    rho.update.details <- rho.update(rho, c.mat, g1.mat, g2.mat, latent,
                                     tau, r, cormat.inv, cormat.logdet, cur.ss.latent,
                                     rho.upper = rho.upper, att.rho, acc.rho, mh.rho)
    rho <- rho.update.details$rho
    cormat.inv <- rho.update.details$cormat.inv
    cormat.logdet <- rho.update.details$cormat.logdet
    cur.ss.latent <- rho.update.details$cur.ss.latent
    att.rho <- rho.update.details$att.rho
    acc.rho <- rho.update.details$acc.rho

    r.update.details <- r.update(ns, nmesh, r, tau, cur.rss, cur.ss.latent, att.r, acc.r, mh.r)
    r <- r.update.details$r
    att.r <- r.update.details$att.r
    acc.r <- r.update.details$acc.r

    Y <- fill.censored.y(Y, cutoff.Y, censored.cases, X.theta, A.latent, tau, r)

    crossprod.Xy <- as.vector(crossprod(X, Y))
    crossprod.Ay <- as.vector(crossprod(A, Y))

    if(iter < (burn / 2)){
      this.update <- mhupdate(acc = acc.rho, att = att.rho, mh = mh.rho)
      acc.rho <- this.update$acc
      att.rho <- this.update$att
      mh.rho <- this.update$mh

      this.update <- mhupdate(acc = acc.r, att = att.r, mh = mh.r)
      acc.r <- this.update$acc
      att.r <- this.update$att
      mh.r <- this.update$mh
    }
  }

    if(iter > burn){
      latent.sum <- latent.sum + latent
      latent2.sum <- latent2.sum + latent^2

      Y.pred <- c(X.pred %*% theta) + as.vector(A.pred %*% latent) + (sqrt(1 - r) * rnorm(np)) / sqrt(tau)

      Y.pred.sum <- Y.pred.sum + Y.pred
      Y.pred2.sum <- Y.pred2.sum + Y.pred^2
    }

    # storage
    keepers.theta[iter, ] <- theta
    keepers.tau[iter] <- tau
    keepers.rho[iter] <- rho
    keepers.r[iter] <- r

  } #end iters

  latent.posmean <- latent.sum / length(return.iters)
  latent.posvar <- latent2.sum / length(return.iters) - latent.posmean^2

  Y.pred.posmean <- Y.pred.sum / length(return.iters)
  Y.pred.posvar <- Y.pred2.sum / length(return.iters) - Y.pred.posmean^2

  tock <- proc.time()[3]

  results <- list(theta = keepers.theta,
                  tau = keepers.tau,
                  rho = keepers.rho,
                  r = keepers.r,
                  latent.posmean = latent.posmean,
                  latent.posvar = latent.posvar,
                  Y.pred.posmean = Y.pred.posmean,
                  Y.pred.posvar = Y.pred.posvar,
                  minutes = (tock - tick) / 60
  )

  return(results)}
