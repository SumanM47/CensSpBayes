#' @name CensSpBayesVS
#' @title Function to run CensSpBayes model on left censored massive spatial data with the option to employ variable selection
#'
#' @description
#' Generate posterior samples for the parameters in the model and predictions of the true process at the prediction locations with pointwise variances for them
#'
#' @usage
#' CensSpBayesVS(Y, S, X, cutoff_Y, S_pred, X_pred, inla_mats, alpha = 2,
#'                theta_init = NULL, latent_init = NULL, tau_init = NULL,
#'                rho_init = 0.5, r_init = NULL,
#'                jitter=0,
#'                mean_theta = 0, sd_theta = 1e2,
#'                tau_a = 0.1, tau_b = 0.1,
#'                rho_upper = NULL,
#'                prior = c("gaussian","hsp"),
#'                vss = c("none","cr","2means","hsp"),
#'                addpar = NULL,
#'                iters = 4000, burn = 2000, thin = 5)
#'
#'
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
#' @param prior character string of length 1, choice of prior. Must be either "gaussian" for the baseline model with Gaussian prior or "hsp" for the model with horseshoe+ prior.
#' @param vss character string of length 1, choice of variable selection strategy. Must be one of "none", "cr", "2means" or "hsp". The option "none" means no post-hoc variable selection will be done, while for "cr', "2means", and "hsp", the variable selection will be done using credible interval method, two-means clustering method, or using horseshoe+ inclusion probability method. Choosing "gaussian" prior and "non" vss reverts to the standard CensSpBayes model.
#' @param addpar additional parameter for variable selection strategy. Unused if vss = "none". For vss = "cr" it is a number between 0 and 1 defining the confidence coefficient. For vss = "2means" it is a positive integer used as a tuning parameter for the two means algorithm. For vss = "hsp" it is a number between 0 and 1 giving the cutoff for posterior inclusion probability for the horseshoe+ model beyond which we are to include a covariate into the model as significant.
#' @param iters number of posterior samples to be drawn
#' @param burn number of posterior samples to be discarded as burn-in period samples
#' @param thin thinning interval. The toatl number of iterations is thin*iters
#'
#'
#' @importFrom Matrix diag
#' @import stats
#' @import spam
#' @import VsusP
#'
#' @return a list of posterior samples for theta, tau, rho and r. Additionally, includes the posterior mean and variances of the latent process and the predicted process, the computation time in minutes, and the indices of covariates deemed significant
#' @export
#'


CensSpBayesVS <- function(Y, S, X, cutoff_Y, S_pred, X_pred, inla_mats, alpha = 2,
                          theta_init = NULL, latent_init = NULL, tau_init = NULL,
                          rho_init = 0.5, r_init = NULL,
                          jitter=0,
                          mean_theta = 0, sd_theta = 1e2,
                          tau_a = 0.1, tau_b = 0.1,
                          rho_upper = NULL,
                          prior = c("gaussian","hsp"),
                          vss = c("none","cr","2means","hsp"),
                          addpar=NULL,
                          iters = 4000, burn = 2000, thin = 5){

  if(is.null(prior)){prior <- "hsp"}
  prior <- tolower(prior)
  if(length(prior)!=1L || !prior %in% c("gaussian","hsp")){stop("prior must be a character of length 1 with values in c('gaussian','hsp')")}

  if(is.null(vss)){vss <- "cr"}
  vss <- tolower(vss)
  if(length(vss)!=1L || !vss %in% c("none","cr","2means","hsp")){stop("vss must be a charcter of length 1 with values in c('none','cr','2means','hsp')")}

  if(prior=="gaussian" & vss == "hsp"){cat("Horsehsoe+ based variable selection can only be done with Horseshoe+ prior. Variable selection strategy switched to 'cr' \n"); vss <- "cr"}

  if(ncol(X)==1){vss <- "none"; cat("Only one variable, no variable selection will be done\n")}

  if(vss=="hsp"){
    if(is.null(addpar)){addpar=0.5}
    if(length(addpar) > 1| addpar >= 1 | addpar <= 0){stop("For vss='hsp', the addpar parameter must be a scalar between 0 and 1")}
    kappacut <- 1-addpar
  }
  if(vss=="cr"){
    if(is.null(addpar)){addpar=0.95}
    if(length(addpar) > 1| addpar >= 1 | addpar <=0){stop("For vss='cr', the addpar parameter must be a scalar between 0 and 1")}
    llcr <- (1-addpar)/2
    ulcr <- 1 - llcr
  }
  if(vss=="2means"){
    if(is.null(addpar)){addpar=20}
    if(length(addpar) > 1 | addpar <=0 | (addpar - ceiling(addpar)) !=0){stop("For vss='2means', the addpar parameter must be a positive integer")}
    l2use <- addpar
  }

  if(prior=="gaussian"){
    out <- CensSpBayes::CensSpBayes(Y, S, X, cutoff_Y, S_pred, X_pred, inla_mats, alpha,
                                    theta_init, latent_init, tau_init,
                                    rho_init, r_init,
                                    jitter,
                                    mean_theta, sd_theta,
                                    tau_a, tau_b,
                                    rho_upper,
                                    iters, burn, thin)
  }

  if(prior=="hsp"){
    out <- CensSpBayes::CensSpBayes2(Y, S, X, cutoff_Y, S_pred, X_pred, inla_mats, alpha,
                                    theta_init, latent_init, tau_init,
                                    rho_init, r_init,
                                    jitter,
                                    mean_theta, sd_theta,
                                    tau_a, tau_b,
                                    rho_upper,
                                    iters, burn, thin)
  }


  if(vss=="cr"){
    out$ind2inc <- which(apply(out$theta[(burn+1):iters,],2,"quantile",probs=llcr)*apply(out$theta[(burn+1):iters,],2,"quantile",probs=ulcr)>0)
  }
  if(vss=="2means"){
    temp <- VsusP::Sequential2MeansBeta(out$theta[(burn+1):iters,],lower=0,upper=1,l=l2use)
    H2use <- as.numeric(names(which.max(table(temp$H.b.i))))
    vsuspinds <- VsusP::S2MVarSelection(out$theta[(burn+1):iters,],H2use)
    out$ind2inc <- vsuspinds
  }
  if(vss=="hsp"){
    out$ind2inc <- which(colMeans(out$kappa[(burn+1):iters,]) < kappacut)
  }
  if(vss=="none"){
    out$ind2inc <- 1:ncol(out$theta)
  }

  bmn <- colMeans(out$theta)
  bmn[-out$ind2inc] <- 0
  Xb <- c(X%*%bmn)
  out$Y_pred_posmean <- Xb + as.numeric(out$inla_mats_used$A_pred%*%c(out$latent_posmean))
  th2use <- out$theta
  th2use[-out$ind2inc] <- 0
  out$theta <- th2use

  return(out)
}


