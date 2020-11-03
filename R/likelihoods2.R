#' Likelihood for M3
#'
#' This model is similar to M2 but allows the inclusions of samples that only
#' have overall expression and genotypes.
#' @export
loglM3 <- function(pars, dat){
  outcs   <- dat$outcs
  gts     <- dat$gts
  tot.exp <- dat$tot.exp
  tot.gts <- dat$tot.gts

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
  la <- dmvnorm(x = outcs-meanmatA, mean = c(0,0), sigma = omega)
  li <- dmvnorm(x = outcs-meanmatI, mean = c(0,0), sigma = omega)
  l <- la*gamma + li*(1-gamma)
  # Samples without Xa/Xi-expression still contain information
  tot.gamma <- exp(besc0+besc*tot.gts)/(1+exp(besc0+besc*tot.gts));
  tot.la <- dnorm(tot.exp, mean = bxa0+bxa*tot.gts +
                                  bxi0+bxi*tot.gts, sd=sqrt(sum(omega)))
  tot.li <- dnorm(tot.exp, mean = bxa0+bxa*tot.gts, sd=sqrt(sum(omega)))
  tot.l <- tot.la*tot.gamma + tot.li*(1-tot.gamma);

  logL <- -sum(log(l))-sum(log(tot.l));
  return(logL);
}

#' @export
loglM3H0_xi <- function(pars, dat){
  pars_H0 <- c(pars[1:3],0,pars[5:length(pars)]);
  ret <- loglM2(pars = pars_H0, dat = dat)
  return(ret)
}
