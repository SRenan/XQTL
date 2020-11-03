# FULL MODEL -------------------
#' @export
loglM2 <- function(pars, dat, finiteonly = TRUE){
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
  meanmatA <- matrix(c(bxa0+bxa*gts, bxi0+bxi*gts), ncol = 2)
  meanmatI <- matrix(c(bxa0+bxa*gts, 0*gts), ncol = 2)
  # Fit
  la <- dmvnorm(x = outcs-meanmatA, mean = c(0,0), sigma = omega, log = T)
  li <- dmvnorm(x = outcs-meanmatI, mean = c(0,0), sigma = omega, log = T)
  l <- exp(la)*gamma + exp(li)*(1-gamma)
  logL <- log(l)
  if(finiteonly){
    ret <- -sum(logL[is.finite(logL)])
  } else{
    ret <- -sum(logL)
  }
  return(ret)
}

# REDUCED MODEL  -------------------
#' @export
loglM2H0_xa <- function(pars, dat, finiteonly = TRUE){
  pars_H0 <- c(pars[1],0,pars[3:length(pars)]);
  ret <- loglM2(pars = pars_H0, dat = dat, finiteonly = finiteonly)
  return(ret)
}
#' @export
loglM2H0_xi <- function(pars, dat, finiteonly = TRUE){
  pars_H0 <- c(pars[1:3],0,pars[5:length(pars)]);
  ret <- loglM2(pars = pars_H0, dat = dat, finiteonly = finiteonly)
  return(ret)
}
