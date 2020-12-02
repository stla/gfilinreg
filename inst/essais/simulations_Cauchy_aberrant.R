library(gfilinreg)
library(heavy)
library(data.table)

nsims <- 500L
MAXLHD <- matrix(NA_real_, nrow = nsims, ncol = 3L)
colnames(MAXLHD) <- c("group1", "group2", "sigma")
FIDlist <- vector("list", length = nsims)

group <- gl(2L, 6L)
set.seed(666L)
for(i in 1L:333L){
  cat(i, " - ")
  # simulated dataset
  dat <- data.frame(
    group = group,
    y = c(rcauchy(6L), 2 + rcauchy(6L))
  )
  # # max-likelihood estimates
  # hfit <- heavyLm(y ~ 0 + group, data = dat, family = Cauchy())
  # MAXLHD[i, ] <- c(hfit[["coefficients"]], sqrt(hfit[["sigma2"]]))
  # # fiducial stuff
  # fidsamples <- gfilinreg(y ~ 0 + group, data = dat, L = 100L, distr = "cauchy")
  # FIDlist[[i]] <- cbind(
  #   parameter = c("group1", "group2", "sigma"), 
  #   as.data.table(gfiSummary(fidsamples))
  # )
}

library(MASS)
X <- model.matrix(~ 0 + group, data = dat)
likelihood <- function(y, beta0, beta1, sigma){
  prod(dcauchy((y-X%*%c(beta0,beta1))/sigma)/sigma)
}
(ML <- MASS::fitdistr(dat$y, likelihood, list(beta0=0, beta1=9500, sigma=5)))
