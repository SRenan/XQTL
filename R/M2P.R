# POISSON -----------------
# FULL MODEL -------------------
#' @export
loglM2P <- function(pars, dat, finiteonly = TRUE){
  # data
  outcs <- dat$expression
  gts   <- dat$genotype
  # parameters
  bxa0     <- exp(pars[1])
  bxa      <- pars[2]
  bxi0     <- exp(pars[3])
  bxi      <- pars[4]
  besc0    <- pars[5]
  besc     <- pars[6]
  varXa    <- exp(pars[7])
  varXi    <- exp(pars[8])
  corXaXi  <- exp(pars[9])/(1+exp(pars[9]))
  covXaXi  <- corXaXi*sqrt(varXa*varXi)
  omega <- matrix(c(varXa, covXaXi, covXaXi, varXi), ncol = 2)
  # Add data
  gamma <- exp(besc0+besc*gts)/(1+exp(besc0+besc*gts))
  #meanmatA <- matrix(c(bxa0+bxa*gts, bxi0+bxi*gts), ncol = 2)
  #meanmatI <- matrix(c(bxa0+bxa*gts, 0*gts), ncol = 2)
  lambdaa <- exp(bxa0+bxa*gts + bxa0+bxi*gts)
  lambdai <- exp(bxa0+bxa*gts + 0*gts)
  # Fit
  la <- dpois(x = outcs[,1], lambda = lambdaa, log = T)
  li <- dpois(x = outcs[,2], lambda = lambdai, log = T)
  # la <- dmvnorm(x = outcs-meanmatA, mean = c(0,0), sigma = omega, log = T)
  # li <- dmvnorm(x = outcs-meanmatI, mean = c(0,0), sigma = omega, log = T)
  l <- exp(la)*gamma + exp(li)*(1-gamma)
  logL <- log(l)
  if(finiteonly){
    ret <- -sum(logL[is.finite(logL)])
  } else{
    ret <- -sum(logL)
  }
  return(ret)
}


# NEGATIVE BINOMIAL ---------------
# FULL MODEL -------------------
#' @importFrom RMKdiscrete dbinegbin
#' @export
loglM2NB <- function(pars, dat, transform = TRUE, finiteonly = TRUE){
  # data
  outcs <- dat$expression
  gts   <- dat$genotype
  if(transform){
    bxs0     <- exp(pars[1])
    bxs      <- pars[2]  # 0 in M0
    bxa0     <- exp(pars[3])
    bxa      <- pars[4]
    bxi0     <- exp(pars[5])
    bxi      <- pars[6]  # 0 in M0
    # p are (0,1]
    p0 <- exp(pars[7])/(1+exp(pars[7]))
    pa <- exp(pars[8])/(1+exp(pars[8]))
    pi <- exp(pars[9])/(1+exp(pars[9]))
  } else{
    bxs0     <- pars[1]
    bxs      <- pars[2]  # 0 in M0
    bxa0     <- pars[3]
    bxa      <- pars[4]
    bxi0     <- pars[5]
    bxi      <- pars[6]  # 0 in M0
    # p are (0,1]
    p0 <- pars[7]
    pa <- pars[8]
    pi <- pars[9]
  }
  ps <- c(p0, pa, pi)
  # Add data
  mu0 <- exp(bxs0+bxs*gts)
  mua <- exp(bxa0+bxa*gts)
  mui <- exp(bxi0+bxi*gts)
  nu0 <- mu0*p0/(1-p0)
  nua <- mua*pa/(1-pa)
  nui <- mui*pi/(1-pi)
  nus <- c(nu0, nua, nui)
  # Fit
  logL <- dbinegbin(y = outcs, nu = nus, p = ps, log = T)

  if(finiteonly){
    ret <- -sum(logL[is.finite(logL)])
  } else{
    ret <- -sum(logL)
  }
  return(ret)
}

# REDUCED MODEL  -------------------
#' @export
loglM2NBH0 <- function(pars, dat, finiteonly = TRUE){
  # Assume no genetic component in the Xi and shared mixtures?
  pars_H0 <- c(pars[1],0, pars[3:5], 0, pars[7:length(pars)]);
  ret <- loglM2NB(pars = pars_H0, dat = dat, finiteonly = finiteonly)
  return(ret)
}
