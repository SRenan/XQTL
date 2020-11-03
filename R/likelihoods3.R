#' Likelihoods for M4
#'
#' M4
#'
#' @importFrom mvtnorm dmvnorm
#' @export
loglM4 <- function(pars, dat, finiteonly = TRUE){
  outcs   <- dat$expression;
  gts     <- dat$genotype;
  tot.exp <- dat$tot.exp;
  tot.gts <- dat$tot.gts;

  ## parameters
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
  la <- dmvnorm(x = outcs-meanmatA, mean = c(0,0), sigma = omega)
  li <- dmvnorm(x = outcs-meanmatI, mean = c(0,0), sigma = omega)
  l <- la*gamma + li*(1-gamma)
  # Samples without Xa/Xi-expression still contain information
  tot.gamma <- exp(besc0+besc*tot.gts)/(1+exp(besc0+besc*tot.gts));
  tot.la <- dnorm(tot.exp, mean = bxa0+bxa*tot.gts +
                                  bxi0+bxi*tot.gts, sd=sqrt(sum(omega)))
  tot.li <- dnorm(tot.exp, mean = bxa0+bxa*tot.gts, sd=sqrt(sum(omega)))
  tot.l <- tot.la*tot.gamma + tot.li*(1-tot.gamma);

  # logL <- -sum(log(l))-sum(log(tot.l));
  logL <- c(log(l), log(tot.l))
  if(finiteonly){
    ret <- -sum(logL[is.finite(logL)])
  } else{
    ret <- -sum(logL)
  }
  return(ret)
}
#' @export
loglM4H0_xi <- function(pars, dat, finiteonly = TRUE) {
  pars_H0 <- c(pars[1:3],0,pars[4:length(pars)]);
  return(loglM4(pars_H0, dat=dat, finiteonly = finiteonly))
}
#' Model M5
#'
#' This is the same model as M4 but with the assumption that bxi = 0. It is used
#' to test besc = 0 under the assumption that there is no effect on Xi. In a
#' broader sense, it should find any QTL that affects expression by modifying
#' escape (chance of escape or amount of escape).
#'
#' @export
loglM5 <- function(pars,dat, finiteonly = TRUE) {
  return(loglM4H0_xi(pars,dat, finiteonly = finiteonly));
}
#' @export
loglM5H0_esc <- function(pars,dat, finiteonly = TRUE) {
  pars.full <- c(pars[1:4],0,pars[5:7]);
  return(loglM5(pars.full, dat, finiteonly = finiteonly))
}

