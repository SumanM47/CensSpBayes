#' @name CensSpBayes3
#' @title Function to run MCMC for left censored massive spatial data with variable selection
#'
#' @description
#' Generate posterior samples for the parameters in the model and predictions of the true process at the prediction locations with pointwise variances for them
#'
#' @usage CensSpBayes3(Y, S, X, cutoff_Y, S_pred, X_pred, inla_mats, alpha = 2,
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

CensSpBayes3  <- function(Y, S, X, cutoff_Y, S_pred, X_pred, inla_mats, alpha = 2,
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

  Y_nc <- Y[-censored_cases]
  Y_c <- Y[censored_cases]

  X_nc <- X[-censored_cases,]
  X_c <- X[censored_cases,]

  A_nc <- A[-censored_cases,]
  A_c <- A[censored_cases,]

  nSc <- length(censored_cases)
  nSnc <- ns - nSc

  # print("CP 1")

  if(is.null(theta_init)){
    theta <- c(solve(crossprod(X)+1e-6*diag(nq)) %*% crossprod(X, Y))
  }else{theta <- theta_init}
  # X_theta <- c(X %*% theta)
  Xc_theta <- c(X_c%*%theta)
  Xnc_theta <- c(X_nc%*%theta)

  if(is.null(latent_init)){
    latent <- rep(0, nmesh)
  }else{latent <- latent_init}
  # Aw <- c(A%*%latent)
  Acw <- as.vector(A_c%*%latent)
  Ancw <- as.vector(A_nc%*%latent)

  if(is.null(tau_init)){
    tau <- 1 / var(Y_nc - Xnc_theta)
  }else{tau <- tau_init}

  if(is.null(r_init)){
    r <- 0.5
  }else{r <- r_init}

  rho <- rho_init
  if(is.null(rho_upper)){rho_upper <- 0.25*max(dist(S))}

  tau_beta <- 0.5

  lambda <- rep(sd_theta,nq)

  L_theta <- 2*ceiling((nq^0.25))
  L_latent <- 2*ceiling((nmesh^0.25))




  Q <- ((rho^2) / (4 * pi))*(((1 / rho)^4) * c_mat + (2 * (1 / rho)^2) * g1_mat + g2_mat)
  Q_ch <- spam::chol(Q)
  Q_logdet <- 2 * sum(log(Matrix::diag(Q_ch)))


  Qlat <- Q%*%latent
  ss_latent <- sum(Qlat*latent)
  ss_nc <- sum((Y_nc - Xnc_theta - Ancw)^2)
  res_c <- Y_c - Xc_theta - Acw
  res_nc <- Y_nc - Xnc_theta - Ancw

  # print("CP 2")

  acc_theta <- acc_latent <- acc_tau <- acc_rho <- acc_r <- acc_taubeta <- att <- 0
  mh_tau <- mh_rho <- mh_r <- mh_taubeta <- 0.5
  # acc_lambda <- rep(0,nq)
  # mh_lambda <- rep(0.1,nq)
  acc_lambda <- 0
  mh_lambda <- 0.1
  mh_theta <- 0.5/L_theta
  mh_latent <- 0.5/L_latent



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

  # print("CP 3")

  for(iter in 1:iters){for(ttt in 1:thin){

    ## Update theta -- HMC

    can_theta <- theta
    can_P <- rnorm(nq)/lambda
    P <- can_P

    ## Compute loglikelihood and gradient

    cur_l_theta <- -0.5*(tau/(1-r))*ss_nc + sum(pnorm(sqrt(tau/(1-r))*res_c,log=T)) - 0.5*sum(((theta - mean_theta)/lambda)^2)
    cur_grad_theta <- (tau/(1-r))*(crossprod(X_nc,res_nc)) - sqrt(tau/(1-r))*(crossprod(X_c,exp(dnorm(sqrt(tau/(1-r))*res_c,log=TRUE)-pnorm(sqrt(tau/(1-r))*res_c,log=TRUE)))) - (theta - mean_theta)/(lambda^2)

    # print(sum(is.na(cur_grad_theta)))

    ## Half-step before the start

    can_P <- can_P - 0.5*mh_theta*cur_grad_theta
    for(iind in 1:L_theta){
      can_theta <- can_theta + mh_theta*can_P*(lambda^2)

      # print(sum(is.na(can_theta)))

      can_Xnc_theta <- c(X_nc%*%can_theta)
      can_Xc_theta <- c(X_c%*%can_theta)
      # can_Acw <- as.vector(A_c%*%can_latent)
      # can_Ancw <- as.vector(A_nc%*%can_latent)

      # can_Qlat <- Q%*%can_latent
      # can_ss_latent <- sum(can_Qlat*can_latent)
      can_res_nc <- Y_nc - can_Xnc_theta - Ancw
      can_ss_nc <- sum(can_res_nc^2)
      can_res_c <- Y_c - can_Xc_theta - Acw

      can_grad_theta <- (tau/(1-r))*(crossprod(X_nc,can_res_nc)) - sqrt(tau/(1-r))*(crossprod(X_c,exp(dnorm(sqrt(tau/(1-r))*can_res_c,log=TRUE)-pnorm(sqrt(tau/(1-r))*can_res_c,log=TRUE)))) - (can_theta - mean_theta)/(lambda^2)

      if(iind!=L_theta){can_P <- can_P - mh_theta*can_grad_theta}
    }

    ## Final half-step

    can_P <- can_P - 0.5*mh_theta*can_grad_theta

    ## Negate momentum

    can_P <- -can_P

    can_l_theta <- -0.5*(tau/(1-r))*can_ss_nc + sum(pnorm(sqrt(tau/(1-r))*can_res_c,log=T)) - 0.5*sum(((can_theta - mean_theta)/lambda)^2)

    # print(c(sum(is.na(can_theta)),sum(is.na(can_P)),can_ss_nc,tau,r,can_l_theta,cur_l_theta))

    a_theta <- can_l_theta - cur_l_theta - 0.5*sum((lambda*can_P)^2) + 0.5*sum((lambda*P)^2)
    # a_theta <- ifelse(is.na(a_theta),-Inf,a_theta)

    if(log(runif(1))<a_theta){
      theta <- can_theta
      Xc_theta <- can_Xc_theta
      Xnc_theta <- can_Xnc_theta
      # Ancw <- can_Ancw
      # Acw <- can_Acw
      # ss_latent <- can_ss_latent
      ss_nc <- can_ss_nc
      res_c <- can_res_c
      res_nc <- can_res_nc
      # Qlat <- can_Qlat
      acc_theta <- acc_theta+1
    }

    ## Update latent -- HMC

    can_latent <- latent
    can_P <- sqrt(r/tau)*as.vector(t(Q_ch)%*%rnorm(nmesh))
    P <- can_P

    ## Compute loglikelihood and gradient

    cur_l_latent <- -0.5*(tau/(1-r))*ss_nc + sum(pnorm(sqrt(tau/(1-r))*res_c,log=T)) - 0.5*(tau/r)*ss_latent
    cur_grad_latent <- (tau/(1-r))*(crossprod(A_nc,res_nc)) - sqrt(tau/(1-r))*(crossprod(A_c,exp(dnorm(sqrt(tau/(1-r))*res_c,log=TRUE)-pnorm(sqrt(tau/(1-r))*res_c,log=TRUE)))) - (tau/r)*Qlat

    ## Half-step before the start

    can_P <- can_P - 0.5*mh_latent*cur_grad_latent
    for(iind in 1:L_latent){
      can_latent <- can_latent + (tau/r)*as.vector(mh_latent*spam::solve(Q,can_P))

      # can_Xnc_theta <- c(X_nc%*%can_theta)
      # can_Xc_theta <- c(X_c%*%can_theta)
      can_Acw <- as.vector(A_c%*%can_latent)
      can_Ancw <- as.vector(A_nc%*%can_latent)

      can_Qlat <- Q%*%can_latent
      can_ss_latent <- sum(can_Qlat*can_latent)
      can_res_nc <- Y_nc - Xnc_theta - can_Ancw
      can_ss_nc <- sum(can_res_nc^2)
      can_res_c <- Y_c - Xc_theta - can_Acw

      can_grad_latent <- (tau/(1-r))*(crossprod(A_nc,can_res_nc)) - sqrt(tau/(1-r))*(crossprod(A_c,exp(dnorm(sqrt(tau/(1-r))*can_res_c,log=TRUE)-pnorm(sqrt(tau/(1-r))*can_res_c,log=TRUE)))) - (tau/r)*can_Qlat

      if(iind!=L_latent){can_P <- can_P - mh_latent*can_grad_latent}
    }

    ## Final half-step

    can_P <- can_P - 0.5*mh_latent*can_grad_latent

    ## Negate momentum

    can_P <- -can_P

    can_l_latent <- -0.5*(tau/(1-r))*can_ss_nc + sum(pnorm(sqrt(tau/(1-r))*can_res_c,log=T)) - 0.5*(tau/r)*can_ss_latent

    a_latent <- can_l_latent - cur_l_latent - 0.5*(tau/r)*sum(spam::backsolve(Q_ch,can_P)^2) + 0.5*(tau/r)*sum(spam::backsolve(Q_ch,P)^2)
    # a_latent <- ifelse(is.na(a_latent),-Inf,a_latent)

    if(log(runif(1))<a_latent){
      latent <- can_latent
      # Xc_theta <- can_Xc_theta
      # Xnc_theta <- can_Xnc_theta
      Ancw <- can_Ancw
      Acw <- can_Acw
      ss_latent <- can_ss_latent
      ss_nc <- can_ss_nc
      res_c <- can_res_c
      res_nc <- can_res_nc
      Qlat <- can_Qlat
      acc_latent <- acc_latent+1
    }



    # cur_l_theta <- -0.5*(tau/(1-r))*ss_nc + sum(pnorm(sqrt(tau/(1-r))*res_c,log=T)) - 0.5*sum(((theta - mean_theta)/lambda)^2)

    # cur_l_thetalatent <- -0.5*(tau/(1-r))*ss_nc + sum(pnorm(sqrt(tau/(1-r))*res_c,log=T)) - 0.5*sum(((theta - mean_theta)/lambda)^2) - 0.5*(tau/r)*ss_latent
    # cur_grad_theta <- (tau/(1-r))*(crossprod(X_nc,res_nc)) - sqrt(tau/(1-r))*(crossprod(X_c,exp(dnorm(sqrt(tau/(1-r))*res_c,log=TRUE)-pnorm(sqrt(tau/(1-r))*res_c,log=TRUE)))) - (theta - mean_theta)/(lambda^2)
    # cur_grad_latent <- (tau/(1-r))*(crossprod(A_nc,res_nc)) - sqrt(tau/(1-r))*(crossprod(A_c,exp(dnorm(sqrt(tau/(1-r))*res_c,log=TRUE)-pnorm(sqrt(tau/(1-r))*res_c,log=TRUE)))) - (tau/r)*Qlat
    #
    #
    # can_theta <- theta + 0.5*mh_thetalatent*cur_grad_theta + sqrt(mh_thetalatent)*rnorm(nq)
    # can_latent <- latent + 0.5*mh_thetalatent*cur_grad_latent + sqrt(mh_thetalatent)*rnorm(nmesh)
    #
    #
    # can_Xnc_theta <- c(X_nc%*%can_theta)
    # can_Xc_theta <- c(X_c%*%can_theta)
    # can_Acw <- as.vector(A_c%*%can_latent)
    # can_Ancw <- as.vector(A_nc%*%can_latent)
    #
    # can_Qlat <- Q%*%can_latent
    # can_ss_latent <- sum(can_Qlat*can_latent)
    # can_ss_nc <- sum((Y_nc - can_Xnc_theta - can_Ancw)^2)
    # can_res_c <- Y_c - can_Xc_theta - can_Acw
    # can_res_nc <- Y_nc - can_Xnc_theta - can_Ancw


    # can_ss_nc <- sum((Y_nc - can_Xnc_theta - Ancw)^2)
    # can_res_c <- Y_c - can_Xc_theta - Acw
    # can_l_theta <- -0.5*(tau/(1-r))*can_ss_nc + sum(pnorm(sqrt(tau/(1-r))*can_res_c,log=T)) - 0.5*sum(((can_theta - mean_theta)/lambda)^2)
    # can_l_thetalatent <- -0.5*(tau/(1-r))*can_ss_nc + sum(pnorm(sqrt(tau/(1-r))*can_res_c,log=T)) - 0.5*sum(((can_theta - mean_theta)/lambda)^2) - 0.5*(tau/r)*can_ss_latent

    # can_grad_theta <- (tau/(1-r))*(crossprod(X_nc,can_res_nc)) - sqrt(tau/(1-r))*(crossprod(X_c,exp(dnorm(sqrt(tau/(1-r))*can_res_c,log=TRUE)-pnorm(sqrt(tau/(1-r))*can_res_c,log=TRUE)))) - (can_theta - mean_theta)/(lambda^2)
    # can_grad_latent <- (tau/(1-r))*(crossprod(A_nc,can_res_nc)) - sqrt(tau/(1-r))*(crossprod(A_c,exp(dnorm(sqrt(tau/(1-r))*can_res_c,log=TRUE)-pnorm(sqrt(tau/(1-r))*can_res_c,log=TRUE)))) - (tau/r)*can_Qlat


    # a_theta <- can_l_theta - cur_l_theta
    # a_thetalatent <- can_l_thetalatent - cur_l_thetalatent - 0.5*(sum((theta - can_theta - 0.5*mh_thetalatent*can_grad_theta)^2)+sum((latent - can_latent - 0.5*mh_thetalatent*can_grad_latent)^2))/mh_thetalatent + 0.5*(sum((can_theta - theta - 0.5*mh_thetalatent*cur_grad_theta)^2)+sum((can_latent - latent - 0.5*mh_thetalatent*cur_grad_latent)^2))/mh_thetalatent

    # if(log(runif(1)) < a_thetalatent){
    #   theta <- can_theta
    #   latent <- can_latent
    #   Xc_theta <- can_Xc_theta
    #   Xnc_theta <- can_Xnc_theta
    #   Ancw <- can_Ancw
    #   Acw <- can_Acw
    #   ss_latent <- can_ss_latent
    #   ss_nc <- can_ss_nc
    #   res_c <- can_res_c
    #   res_nc <- can_res_nc
    #   Qlat <- can_Qlat
    #   acc_thetalatent <- acc_thetalatent + 1
    # }

    # print("CP 4")


    ## Update latent -- MH (MALA?) --- Not Needed

    # cur_l_latent <- -0.5*(tau/(1-r))*ss_nc + sum(pnorm(sqrt(tau/(1-r))*res_c,log=T)) -0.5*(tau/r)*ss_latent
    #
    #
    # can_latent <- latent + mh_latent*rnorm(nmesh)
    # can_Acw <- as.vector(A_c%*%can_latent)
    # can_Ancw <- as.vector(A_nc%*%can_latent)
    #
    # can_ss_latent <- sum((Q%*%can_latent)*can_latent)
    # can_ss_nc <- sum((Y_nc - Xnc_theta - can_Ancw)^2)
    # can_res_c <- Y_c - Xc_theta - can_Acw
    #
    # can_l_latent <- -0.5*(tau/(1-r))*can_ss_nc + sum(pnorm(sqrt(tau/(1-r))*can_res_c,log=T)) -0.5*(tau/r)*can_ss_latent
    #
    #
    # a_latent <- can_l_latent - cur_l_latent
    #
    # if(log(runif(1)) < a_latent){
    #   latent <- can_latent
    #   Ancw <- can_Ancw
    #   Acw <- can_Acw
    #   ss_latent <- can_ss_latent
    #   ss_nc <- can_ss_nc
    #   res_c <- can_res_c
    #   acc_latent <- acc_latent + 1
    # }


    ## Update tau -- MH

    cur_l_tau <- (0.5*nSnc + 0.5*nmesh + tau_a -1)*log(tau) - ((0.5*ss_nc/(1-r)) + (0.5*ss_latent/r) + tau_b)*tau + sum(pnorm(sqrt(tau/(1-r))*res_c,log=T))

    can_tau <- exp(log(tau) + mh_tau*rnorm(1))

    can_l_tau <- (0.5*nSnc + 0.5*nmesh + tau_a -1)*log(can_tau) - ((0.5*ss_nc/(1-r)) + (0.5*ss_latent/r) + tau_b)*can_tau + sum(pnorm(sqrt(can_tau/(1-r))*res_c,log=T))

    a_tau <- can_l_tau - cur_l_tau + log(can_tau) - log(tau)

    if(log(runif(1)) < a_tau){
      tau <- can_tau
      acc_tau <- acc_tau + 1
    }

    # print("CP 5")

    ## Update rho -- MH

    cur_l_rho <- -0.5*(tau/r)*ss_latent + 0.5*Q_logdet

    rhotr <- rho/rho_upper
    rho_star <- log(rhotr/(1-rhotr))
    can_rhostar <- rho_star + mh_rho*rnorm(1)
    can_rhotr <- 1/(1+exp(-can_rhostar))
    can_rho <- rho_upper*can_rhotr

    can_Q <- ((can_rho^2)/(4*pi))*(((1/can_rho)^4)*c_mat + (2*(1/can_rho)^2) * g1_mat + g2_mat)
    can_Q_ch <- spam::chol(can_Q)
    can_Q_logdet <- 2*sum(log(Matrix::diag(can_Q_ch)))

    can_Qlat <- can_Q%*%latent
    can_ss_latent <- sum(can_Qlat*latent)

    can_l_rho <- -0.5*(tau/r)*can_ss_latent + 0.5*can_Q_logdet

    # print(c(can_ss_latent,can_Q_logdet,can_rho,tau,r))

    a_rho <- can_l_rho - cur_l_rho + log(can_rho) + log(rho_upper - can_rho) - log(rho) - log(rho_upper - rho)

    if(log(runif(1)) < a_rho){
      rho <- can_rho
      Q <- can_Q
      Q_ch <- can_Q_ch
      Q_logdet <- can_Q_logdet
      Qlat <- can_Qlat
      ss_latent <- can_ss_latent
      acc_rho <- acc_rho + 1
    }

    # print("CP 6")


    ## Update r -- MH

    cur_l_r <- -0.5*nSnc*log(1-r) -0.5*nmesh*log(r) -0.5*tau*((ss_nc/(1-r)) + (ss_latent/r)) + sum(pnorm(sqrt(tau/(1-r))*res_c,log=T))

    rtr <- (r - 1e-4)/(0.9999 - 1e-4)
    r_star <- log(rtr/(1-rtr))
    can_rstar <- r_star + mh_r*rnorm(1)
    can_rtr <- 1/(1+exp(-can_rstar))
    can_r <- (0.9999 - 1e-4)*can_rtr + 1e-4

    can_l_r <- -0.5*nSnc*log(1-can_r) -0.5*nmesh*log(can_r) -0.5*tau*((ss_nc/(1-can_r)) + (ss_latent/can_r)) + sum(pnorm(sqrt(tau/(1-can_r))*res_c,log=T))

    a_r <- can_l_r - cur_l_r + log(can_r - 1e-4) + log(0.9999 - can_r) - log(r - 1e-4) - log(0.9999 - r)

    if(log(runif(1)) < a_r){
      r <- can_r
      acc_r <- acc_r + 1
    }

    # print("CP 7")

    ## Update lambda -- MH

    # for(ii in 1:nq){
    #   cur_l_lambda <- -0.5*(((theta[ii] - mean_theta)/lambda[ii])^2) - log(lambda[ii]) + log(abs(log(lambda[ii]/tau_beta))) - log(abs(((lambda[ii]/tau_beta)^2)-1))
    #
    #   # print(cur_l_lambda)
    #
    #   can_lambda <- exp(log(lambda[ii]) + mh_lambda[ii]*rnorm(1))
    #   # print(can_lambda)
    #
    #   can_l_lambda <- -0.5*(((theta[ii] - mean_theta)/can_lambda)^2) - log(can_lambda) + log(abs(log(can_lambda/tau_beta))) - log(abs(((can_lambda/tau_beta)^2)-1))
    #   # print(can_l_lambda)
    #
    #
    #   ratio_lambda <- can_l_lambda - cur_l_lambda + log(can_lambda) - log(lambda[ii])
    #
    #   if(log(runif(1)) < ratio_lambda){
    #     lambda[ii] <- can_lambda
    #     acc_lambda[ii] <- acc_lambda[ii] + 1
    #   }
    # }

    ## Update log lambda -- MALA

    # cur_l_lambda <- -0.5*sum(((theta - mean_theta)/lambda)^2) - sum(log(lambda)) + sum(log(abs(log(lambda/tau_beta)))) - sum(log(abs(((lambda/tau_beta)^2)-1)))
    # cur_grad_lambda <- (((theta - mean_theta)/lambda)^2)/lambda - 1/lambda + 1/(lambda*log(lambda/tau_beta)) - (2*lambda)/((lambda^2) - (tau_beta^2))

    cur_l_lambda <- -0.5*sum(((theta - mean_theta)/lambda)^2) + sum(log(abs(log(lambda/tau_beta)))) - sum(log(abs(((lambda/tau_beta)^2)-1)))
    cur_grad_lambda <- (((theta - mean_theta)/lambda)^2)/lambda + 1/(lambda*log(lambda/tau_beta)) - (2*lambda)/((lambda^2) - (tau_beta^2))

    can_lambda <- exp(log(lambda) + 0.5*mh_lambda*cur_grad_lambda + sqrt(mh_lambda)*rnorm(nq))

    can_l_lambda <- -0.5*sum(((theta - mean_theta)/can_lambda)^2) + sum(log(abs(log(can_lambda/tau_beta)))) - sum(log(abs(((can_lambda/tau_beta)^2)-1)))
    can_grad_lambda <- (((theta - mean_theta)/can_lambda)^2)/lambda + 1/(can_lambda*log(can_lambda/tau_beta)) - (2*can_lambda)/((can_lambda^2) - (tau_beta^2))

    cur_q <- -0.5*sum((can_lambda - lambda - 0.5*mh_lambda*cur_grad_lambda)^2)/mh_lambda
    can_q <- -0.5*sum((lambda - can_lambda - 0.5*mh_lambda*can_grad_lambda)^2)/mh_lambda

    a_lambda <- can_l_lambda + can_q - cur_l_lambda - cur_q
    if(log(runif(1)) < a_lambda){
      lambda <- can_lambda
      acc_lambda <- acc_lambda + 1
    }



    # print("CP 8")
    ## Update tau_beta -- MH

    cur_l_taubeta <- -nq*log(tau_beta) + sum(log(abs(log(lambda/tau_beta)))) - sum(log(abs(((lambda/tau_beta)^2)-1)))

    tau_betastar <- (tau_beta-1e-4)/(0.9999-1e-4)
    can_tau_betastar <- log(tau_betastar/(1-tau_betastar)) + mh_taubeta*rnorm(1)
    can_tau_beta <- 1e-4 + (0.9999-1e-4)/(1+exp(-can_tau_betastar))

    can_l_taubeta <- -nq*log(can_tau_beta) + sum(log(abs(log(lambda/can_tau_beta)))) - sum(log(abs(((lambda/can_tau_beta)^2)-1)))

    ratio_taubeta <- can_l_taubeta - cur_l_taubeta + nq*(log(can_tau_beta) + log(1-can_tau_beta) - log(tau_beta) - log(1-tau_beta))

    if(log(runif(1)) < ratio_taubeta){
      tau_beta <- can_tau_beta
      acc_taubeta <- acc_taubeta + 1
    }

    # print("CP 9")

    att <- att + 1

    if(iter < (burn / 2)){
      if(att > 50){
        acc_rate_theta <- acc_theta/att
        acc_rate_latent <- acc_latent/att
        # acc_rate_thetalatent <- acc_thetalatent/att
        acc_rate_tau <- acc_tau/att
        acc_rate_rho <- acc_rho/att
        acc_rate_r <- acc_r/att
        acc_rate_lambda <- acc_lambda/att
        acc_rate_taubeta <- acc_taubeta/att
        mh_theta <- ifelse(acc_rate_theta < 0.5,mh_theta*0.8,ifelse(acc_rate_theta > 0.7,mh_theta*1.2,mh_theta))
        mh_latent <- ifelse(acc_rate_latent < 0.5,mh_latent*0.8,ifelse(acc_rate_latent > 0.7,mh_latent*1.2,mh_latent))
        # mh_thetalatent <- ifelse(acc_rate_thetalatent < 0.3,mh_thetalatent*0.8,ifelse(acc_rate_thetalatent > 0.5,mh_thetalatent*1.2,mh_thetalatent))
        mh_tau <- ifelse(acc_rate_tau < 0.3,mh_tau*0.8,ifelse(acc_rate_tau > 0.5,mh_tau*1.2,mh_tau))
        mh_rho <- ifelse(acc_rate_rho < 0.3,mh_rho*0.8,ifelse(acc_rate_rho > 0.5,mh_rho*1.2,mh_rho))
        mh_r <- ifelse(acc_rate_r < 0.3,mh_r*0.8,ifelse(acc_rate_r > 0.5,mh_r*1.2,mh_r))
        mh_lambda <- ifelse(acc_rate_lambda < 0.3,mh_lambda*0.8,ifelse(acc_rate_lambda > 0.5,mh_lambda*1.2,mh_lambda))
        mh_taubeta <- ifelse(acc_rate_taubeta < 0.3,mh_taubeta*0.8,ifelse(acc_rate_taubeta > 0.5,mh_taubeta*1.2,mh_taubeta))


        acc_theta <- acc_latent <- acc_tau <- acc_rho <- acc_r <- acc_taubeta <- 0
        # acc_thetalatent <- acc_tau <- acc_rho <- acc_r <- acc_taubeta <- 0
        # acc_lambda <- rep(0,nq)
        acc_lambda <- 0
        att <- 0

      }
    }
  }

    # print("CP 10")

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

    print(iter)

    # print("CP 11")

  } #end iters

  latent_posmean <- latent_sum / length(return_iters)
  latent_posvar <- latent2_sum / length(return_iters) - latent_posmean^2

  Y_pred_posmean <- Y_pred_sum / length(return_iters)
  Y_pred_posvar <- Y_pred2_sum / length(return_iters) - Y_pred_posmean^2

  # propinc <- colMeans(1/(1+(keepers_lambda*keepers_taubeta)^2))
  # incind <- which(propinc < 0.5)
  kapp <- 1/(1+(keepers_lambda*keepers_taubeta)^2)

  tock <- proc.time()[3]

  # print("CP 12")

  results <- list(theta = keepers_theta,
                  tau = keepers_tau,
                  rho = keepers_rho,
                  r = keepers_r,
                  lambda = keepers_lambda,
                  tau_beta = keepers_taubeta,
                  latent_posmean = latent_posmean,
                  latent_posvar = latent_posvar,
                  Y_pred_posmean = Y_pred_posmean,
                  Y_pred_posvar = Y_pred_posvar,
                  kappa = kapp,
                  inla_mats_used = inla_mats,
                  minutes = (tock - tick) / 60
  )

  return(results)}
