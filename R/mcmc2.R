#' @name CensSpBayes2
#' @title Function to run MCMC for left censored massive spatial data with variable selection
#'
#' @description
#' Generate posterior samples for the parameters in the model and predictions of the true process at the prediction locations with pointwise variances for them
#'
#' @usage CensSpBayes2(Y, S, X, cutoff_Y, S_pred, X_pred, inla_mats, alpha = 2,
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

CensSpBayes2 <- function(Y, S, X, cutoff_Y, S_pred, X_pred, inla_mats, alpha = 2,
                        theta_init = NULL, latent_init = NULL, tau_init = NULL,
                        rho_init = 0.5, r_init = NULL,
                        jitter=0,
                        # priors
                        mean_theta = 0, sd_theta = 1,
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

  tau_beta <- 0.5

  lambda <- rep(sd_theta,nq)




  cormat_inv <- ((rho^2) / (4 * pi))*(((1 / rho)^4) * c_mat + (2 * (1 / rho)^2) * g1_mat + g2_mat)
  cormat_inv_c <- spam::chol(cormat_inv)
  cormat_logdet <- -2 * sum(log(Matrix::diag(cormat_inv_c)))


  crossprod_X <- crossprod(X)
  crossprod_A <- crossprod(A)
  crossprod_A_X <- crossprod(A, X)
  crossprod_X_A <- crossprod(X, A)
  crossprod_Xy <- as.vector(crossprod(X, Y))
  crossprod_Ay <- as.vector(crossprod(A, Y))


  acc_rho <- acc_r <- acc_taubeta <- att <- 0
  mh_rho <- mh_r <- mh_taubeta <- 1
  acc_lambda <- rep(0,nq)
  mh_lambda <- rep(1,nq)

  latent_sum <- rep(0, nmesh)
  latent2_sum <- rep(0, nmesh)

  Y_pred_sum <- rep(0, np)
  Y_pred2_sum <- rep(0, np)

  keepers_theta <- matrix(NA, nrow = iters, ncol = nq)
  keepers_tau <- rep(NA, iters)
  keepers_rho <- rep(NA, iters)
  keepers_r <- rep(NA, iters)
  keepers_lambda <- matrix(NA,iters,nq)
  keepers_taubeta <- rep(NA,iters)

  if(jitter > 0){r <- jitter(r); rho <- jitter(rho); tau <- jitter(tau); theta <- jitter(theta)}

  return_iters <- (burn + 1):iters

  for(iter in 1:iters){for(ttt in 1:thin){

    ## Update theta and latent together -- Gibbs

    latent_cov_inv_11 <- tau*crossprod_X / (1 - r) + diag(1/lambda^2)
    latent_cov_inv_12 <- tau*crossprod_X_A / (1 - r)
    latent_cov_inv_21 <- tau*crossprod_A_X / (1 - r)
    latent_cov_inv_22 <- tau*(crossprod_A / (1 - r) + cormat_inv / r)

    latent_mean_part1 <- tau*crossprod_Xy / (1 - r) + mean_theta/(lambda^2)
    latent_mean_part2 <- tau*crossprod_Ay / (1 - r)

    latent_cov_inv <- rbind(cbind(latent_cov_inv_11, latent_cov_inv_12),
                            cbind(latent_cov_inv_21, latent_cov_inv_22))
    latent_mean_part <- c(latent_mean_part1, latent_mean_part2)

    chol_latent_cov_inv <- spam::chol(latent_cov_inv)
    omega <- spam::forwardsolve(t(chol_latent_cov_inv), latent_mean_part)
    mm <- spam::backsolve(chol_latent_cov_inv, omega)
    zz <- rnorm(nq + nmesh)
    vv <- spam::backsolve(chol_latent_cov_inv, zz)
    theta_latent <- mm + vv


    theta <- theta_latent[1:nq]
    latent <- theta_latent[-(1:nq)]
    X_theta <- as.vector(X %*% theta)
    A_latent <- as.vector(A %*% latent)
    cur_rss <- sum((Y - X_theta - A_latent)^2)
    cur_ss_theta <- sum((theta/lambda)^2)
    cur_ss_latent <- sum((cormat_inv %*% latent) * latent)


    ## Update tau -- Gibbs

    tau <- rgamma(1, shape = tau_a + (ns + nmesh) / 2,
                  rate = tau_b + (cur_rss / (1 - r) + cur_ss_latent / r) / 2)


    ## Update rho -- MH

    rho_star <- log((rho/rho_upper)/(1-(rho/rho_upper)))
    can_rho_star <- rnorm(1, rho_star, mh_rho)
    can_rho <- (1/(1+exp(-can_rho_star)))*rho_upper

    can_cormat_inv <- ((can_rho^2) / (4*pi))*(((1/can_rho)^4)*c_mat + (2*(1/can_rho)^2)*g1_mat + g2_mat)
    can_cormat_inv_c <- spam::chol(can_cormat_inv)
    can_cormat_logdet <- -2*sum(log(Matrix::diag(cormat_inv_c)))

    can_ss_latent <- sum((can_cormat_inv %*% latent) * latent)

    ratio_rho <- -0.5 * tau * (can_ss_latent - cur_ss_latent) / r - 0.5 *
      (can_cormat_logdet - cormat_logdet) +
      log(can_rho - 0) + log(rho_upper - can_rho) -
      log(rho - 0) - log(rho_upper - rho)

    if(log(runif(1)) < ratio_rho){
      rho <- can_rho
      cormat_inv <- can_cormat_inv
      cormat_logdet <- can_cormat_logdet
      acc_rho <- acc_rho + 1
      cur_ss_latent <- can_ss_latent}


    ## Update r -- MH

    rtemp <- (r-1e-4)/0.9999
    r_star <- log(rtemp/(1-rtemp))
    can_r_star <- rnorm(1, r_star, mh_r)
    can_rtemp <- 1/(1+exp(-can_r_star))
    can_r <- can_rtemp*(0.9999 - 1e-4) + 1e-4

    ratio1 <- -0.5 * nmesh * (log(can_r) - log(r)) - 0.5 * tau * cur_ss_latent * (1 / can_r - 1 / r)
    ratio2 <- -0.5 * ns * (log(1 - can_r) - log(1 - r)) - 0.5 * tau * cur_rss * (1 / (1 - can_r) - 1 / (1 - r))
    ratio_r <- ratio1 + ratio2 + log(can_r - 1e-4) + log(0.9999 - can_r) - log(r - 1e-4) - log(0.9999 - r)

    if(log(runif(1)) < ratio_r){
      r <- can_r
      acc_r <- acc_r + 1}

    ## Update lambda -- MH

    for(ii in 1:nq){
      cur_l_lambda <- -0.5*(((theta[ii] - mean_theta)/lambda[ii])^2) - log(lambda[ii]) + log(abs(log(lambda[ii]/tau_beta))) - log(abs(((lambda[ii]/tau_beta)^2)-1))

      can_lambda <- exp(log(lambda[ii]) + mh_lambda[ii]*rnorm(1))

      can_l_lambda <- -0.5*(((theta[ii] - mean_theta)/can_lambda[ii])^2) - log(can_lambda[ii]) + log(abs(log(can_lambda[ii]/tau_beta))) - log(abs(((can_lambda[ii]/tau_beta)^2)-1))

      ratio_lambda <- can_l_lambda - cur_l_lambda + log(can_lambda[ii]) - log(lambda[ii])

      if(log(runif(1)) < ratio_lambda){
        lambda[ii] <- can_lambda[ii]
        acc_lambda[ii] <- acc_lambda[ii] + 1
      }
    }

    ## Update tau_beta -- MH

    cur_l_taubeta <- -nq*log(tau_beta) + sum(log(abs(log(lambda/tau_beta)))) - sum(log(abs(((lambda/tau_beta)^2)-1)))

    can_tau_betastar <- log(tau_beta/(1-tau_beta)) + mh_taubeta*rnorm(1)
    can_tau_beta <- 1/(1+exp(-can_tau_betastar))

    can_l_taubeta <- -nq*log(can_tau_beta) + sum(log(abs(log(lambda/can_tau_beta)))) - sum(log(abs(((lambda/can_tau_beta)^2)-1)))

    ratio_taubeta <- can_l_taubeta - cur_l_taubeta + nq*(log(can_tau_beta) + log(1-can_tau_beta) - log(tau_beta) - log(1-tau_beta))

    if(log(runif(1)) < ratio_taubeta){
      tau_beta <- can_tau_beta
      acc_taubeta <- acc_taubeta + 1
    }


    ## Fill Y

    total_mean <- X_theta + A_latent
    excess <- cutoff_Y - total_mean
    total_mean <- total_mean[censored_cases]
    excess <- excess[censored_cases]

    runif_probs <- pnorm(sqrt(tau / (1 - r)) * excess) * runif(length(censored_cases))
    Y[censored_cases] <- total_mean + sqrt((1 - r) / tau) * qnorm(runif_probs)

    crossprod_Xy <- as.vector(crossprod(X, Y))
    crossprod_Ay <- as.vector(crossprod(A, Y))

    att <- att + 1

    if(iter < (burn / 2)){
      if(att > 50){
        acc_rate_rho <- acc_rho/att
        acc_rate_r <- acc_r/att
        acc_rate_lambda <- acc_lambda/att
        acc_rate_taubeta <- acc_taubeta/att
        mh_rho <- ifelse(acc_rate_rho < 0.3,mh_rho*0.8,ifelse(acc_rate_rho > 0.5,mh_rho*1.2,mh_rho))
        mh_r <- ifelse(acc_rate_r < 0.3,mh_r*0.8,ifelse(acc_rate_r > 0.5,mh_r*1.2,mh_r))
        mh_lambda <- ifelse(acc_rate_lambda < 0.3,mh_lambda*0.8,ifelse(acc_rate_lambda > 0.5,mh_lambda*1.2,mh_lambda))
        mh_taubeta <- ifelse(acc_rate_taubeta < 0.3,mh_taubeta*0.8,ifelse(acc_rate_taubeta > 0.5,mh_taubeta*1.2,mh_taubeta))


        acc_rho <- acc_r <- acc_taubeta <- 0
        acc_lambda <- rep(0,nq)
        att <- 0

      }
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
    keepers_lambda[iter,] <- lambda
    keepers_taubeta[iter] <- tau_beta

  } #end iters

  latent_posmean <- latent_sum / length(return_iters)
  latent_posvar <- latent2_sum / length(return_iters) - latent_posmean^2

  Y_pred_posmean <- Y_pred_sum / length(return_iters)
  Y_pred_posvar <- Y_pred2_sum / length(return_iters) - Y_pred_posmean^2

  propinc <- colMeans(1/(1+(keepers_lambda*keepers_taubeta)^2))
  incind <- which(propinc < 0.5)

  tock <- proc.time()[3]

  results <- list(theta = keepers_theta,
                  tau = keepers_tau,
                  rho = keepers_rho,
                  r = keepers_r,
                  latent_posmean = latent_posmean,
                  latent_posvar = latent_posvar,
                  Y_pred_posmean = Y_pred_posmean,
                  Y_pred_posvar = Y_pred_posvar,
                  inclusion_index = incind,
                  minutes = (tock - tick) / 60
  )

  return(results)}
