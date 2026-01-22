#' @name CensSpBayes
#' @title Function to run MCMC for left censored massive spatial data
#'
#' @description
#' Generate posterior samples for the parameters in the model and predictions of the true process at the prediction locations with pointwise variances for them
#'
#' @usage CensSpBayes(Y, S, X, cutoff_Y, S_pred, X_pred, inla_mats, alpha = 2,
#'                theta_init = NULL, latent_init = NULL, tau_init = NULL,
#'                rho_init = 0.5, r_init = NULL,
#'                jitter=0,
#'                mean_theta = 0, sd_theta = 1e2,
#'                tau_a = 0.1, tau_b = 0.1,
#'                rho_upper = NULL,
#'                iters = 4000, burn = 2000, thin = 5)
#'
#' @param Y vector of observations
#' @param S matrix of observation locations with each row corresponding to one location
#' @param X design matrix for covariates, includes intercept
#' @param cutoff_Y vector of cutoff values below which data are censored
#' @param S_pred matrix of prediction locations, each row corresponds to a location
#' @param X_pred matrix of covariate values at the prediction locations
#' @param inla_mats list of matrices needed for the SPDE approximation. The entries are: c.mat(mass matrix), g1.mat(stiffness matrix), g2.mat(a combination of mass and stiffness matrices), A(index matrix for observation locations), A.pred(index matrix for prediction locations)
#' @param alpha the order of the SPDE, defaults to 2, currently the only accepted value
#' @param theta_init initial values for theta, the coefficients for X
#' @param latent_init initial values for latent parameters
#' @param tau_init initial value for the precision parameter
#' @param rho_init initial value for the spatial range parameter
#' @param r_init initial value for the r parameter, the proportion of partial sill to sill
#' @param jitter numeric variable to indicate whether the initial values are to be used exactly (jitter=0) or to start with a slightly different value (jitter !=0)
#' @param mean_theta hyperparameter corresponding to the mean of the normal prior for theta
#' @param sd_theta hyperparameter corresponding to the sd of the normal prior for theta
#' @param tau_a hyperparameter corresponding to the shape of the inverse gamma prior for tau
#' @param tau_b hyperparameter corresponding to the rate of the inverse gamma prior for tau
#' @param rho_upper hyperparameter corresponding to the upper limit of the uniform prior for rho. Recommend use no more than 0.25*maximum distance of the domain
#' @param iters number of posterior samples to be drawn
#' @param burn number of posterior samples to be discarded as burn-in period samples
#' @param thin thinning interval. The toatl number of iterations is thin*iters
#'
#' @importFrom Matrix diag
#' @import stats
#' @import spam
#'
#' @return a list of posterior samples for theta, tau, rho and r. Additionally, includes the posterior mean and variances of the latent process and the predicted process, and the computation time in minutes
#' @export

CensSpBayes <- function(Y, S, X, cutoff_Y, S_pred, X_pred, inla_mats, alpha = 2,
                    theta_init = NULL, latent_init = NULL, tau_init = NULL,
                    rho_init = 0.5, r_init = NULL,
		                jitter=0,
                    # priors
                    mean_theta = 0, sd_theta = 1e2,
                    tau_a = 0.1, tau_b = 0.1,
                    rho_upper = NULL,
                    # mcmc settings
                    iters = 4000, burn = 2000, thin = 5){

  tick <- proc.time()[3]


  c_mat <- inla_mats$c_mat
  g1_mat <- inla_mats$g1_mat
  g2_mat <- inla_mats$g2_mat
  A <- inla_mats$A
  A_pred <- inla_mats$A_pred

  X <- as.matrix(X)
  X_pred <- as.matrix(X_pred)

  ns <- length(Y)
  nmesh <- ncol(A)
  nq <- ncol(X)
  np <- nrow(S_pred)

  censored_cases <- which(Y <= cutoff_Y)
  Y[censored_cases] <- cutoff_Y[censored_cases]

  if(is.null(theta_init)){
    theta <- c(solve(crossprod(X)+1e-6*diag(nq)) %*% crossprod(X, Y))
  }else{theta <- theta_init}
  X_theta <- c(X %*% theta)

  if(is.null(latent_init)){
    latent <- rep(0, nmesh)
  }else{latent <- latent_init}

  if(is.null(tau_init)){
    tau <- 1 / var(Y - X_theta)
  }else{tau <- tau_init}

  if(is.null(r_init)){
    r <- 0.5
  }else{r <- r_init}

  rho <- rho_init
  if(is.null(rho_upper)){rho_upper <- 0.25*max(dist(S))}


  cormat_details <- cormat.inv.update.inla(rho, c_mat, g1_mat, g2_mat, alpha = 2)

  cormat_inv <- cormat_details$cormat.inv
  cormat_logdet <- cormat_details$cormat.logdet

  crossprod_X <- crossprod(X)
  crossprod_A <- crossprod(A)
  crossprod_A_X <- crossprod(A, X)
  crossprod_X_A <- crossprod(X, A)
  crossprod_Xy <- as.vector(crossprod(X, Y))
  crossprod_Ay <- as.vector(crossprod(A, Y))


  acc_rho <- att_rho <- mh_rho <- 1
  acc_r <- att_r <- mh_r <- 1

  latent_sum <- rep(0, nmesh)
  latent2_sum <- rep(0, nmesh)

  Y_pred_sum <- rep(0, np)
  Y_pred2_sum <- rep(0, np)

  keepers_theta <- matrix(NA, nrow = iters, ncol = nq)
  keepers_tau <- rep(NA, iters)
  keepers_rho <- rep(NA, iters)
  keepers_r <- rep(NA, iters)

  if(jitter > 0){r <- jitter(r); rho <- jitter(rho); tau <- jitter(tau); theta <- jitter(theta)}

  return_iters <- (burn + 1):iters

  for(iter in 1:iters){for(ttt in 1:thin){
    theta_latent <- theta.latent.update(nq, nmesh, crossprod_X, crossprod_A,
                                        crossprod_A_X, crossprod_X_A,
                                        crossprod_Xy, crossprod_Ay,
                                        tau, cormat_inv, r, sd_theta)
    theta <- theta_latent[1:nq]
    latent <- theta_latent[-(1:nq)]
    X_theta <- as.vector(X %*% theta)
    A_latent <- as.vector(A %*% latent)
    cur_rss <- sum((Y - X_theta - A_latent)^2)
    cur_ss_theta <- sum(theta^2) / sd_theta^2
    cur_ss_latent <- sum((cormat_inv %*% latent) * latent)

    tau <- tau.update(ns, nq, nmesh, cur_rss, cur_ss_theta, cur_ss_latent, r, tau_a, tau_b)

    rho_update_details <- rho.update(rho, c_mat, g1_mat, g2_mat, latent,
                                     tau, r, cormat_inv, cormat_logdet, cur_ss_latent,
                                     rho.upper = rho_upper, att_rho, acc_rho, mh_rho)
    rho <- rho_update_details$rho
    cormat_inv <- rho_update_details$cormat.inv
    cormat_logdet <- rho_update_details$cormat.logdet
    cur_ss_latent <- rho_update_details$cur.ss.latent
    att_rho <- rho_update_details$att.rho
    acc_rho <- rho_update_details$acc.rho

    r_update_details <- r.update(ns, nmesh, r, tau, cur_rss, cur_ss_latent, att_r, acc_r, mh_r)
    r <- r_update_details$r
    att_r <- r_update_details$att.r
    acc_r <- r_update_details$acc.r

    Y <- fill.censored.y(Y, cutoff_Y, censored_cases, X_theta, A_latent, tau, r)

    crossprod_Xy <- as.vector(crossprod(X, Y))
    crossprod_Ay <- as.vector(crossprod(A, Y))

    if(iter < (burn / 2)){
      this_update <- mhupdate(acc = acc_rho, att = att_rho, mh = mh_rho)
      acc_rho <- this_update$acc
      att_rho <- this_update$att
      mh_rho <- this_update$mh

      this_update <- mhupdate(acc = acc_r, att = att_r, mh = mh_r)
      acc_r <- this_update$acc
      att_r <- this_update$att
      mh_r <- this_update$mh
    }
  }

    if(iter > burn){
      latent_sum <- latent_sum + latent
      latent2_sum <- latent2_sum + latent^2

      Y_pred <- c(X_pred %*% theta) + as.vector(A_pred %*% latent) + (sqrt(1 - r) * rnorm(np)) / sqrt(tau)

      Y_pred_sum <- Y_pred_sum + Y_pred
      Y_pred2_sum <- Y_pred2_sum + Y_pred^2
    }

    # storage
    keepers_theta[iter, ] <- theta
    keepers_tau[iter] <- tau
    keepers_rho[iter] <- rho
    keepers_r[iter] <- r

  } #end iters

  latent_posmean <- latent_sum / length(return_iters)
  latent_posvar <- latent2_sum / length(return_iters) - latent_posmean^2

  Y_pred_posmean <- Y_pred_sum / length(return_iters)
  Y_pred_posvar <- Y_pred2_sum / length(return_iters) - Y_pred_posmean^2

  tock <- proc.time()[3]

  results <- list(theta = keepers_theta,
                  tau = keepers_tau,
                  rho = keepers_rho,
                  r = keepers_r,
                  latent_posmean = latent_posmean,
                  latent_posvar = latent_posvar,
                  Y_pred_posmean = Y_pred_posmean,
                  Y_pred_posvar = Y_pred_posvar,
                  inla_mats_used = inla_mats,
                  minutes = (tock - tick) / 60
  )

  return(results)}
