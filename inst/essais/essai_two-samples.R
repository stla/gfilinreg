library(gfilinreg)

dat <- data.frame(
  group = gl(2, 3),
  y = c(rnorm(3), 2 + rnorm(3))
)
system.time(
  gfi <- gfilinreg(y ~ 0 + group, data = dat, L = 80L, distr = "cauchy")
)

################################################################################
set.seed(666)
dat <- data.frame(
  group = gl(2, 20),
  y = c(rnorm(20), 2 + rnorm(20))
)

gfi <- gfilinreg(y ~ 0 + group, data = dat, L = 100L, distr = "normal")
gfiSummary(gfi)
confint(lm(y ~ 0 + group, data = dat))

gfi <- gfilinreg(y ~ 0 + group, data = dat, L = 100L, distr = "cauchy")
gfiSummary(gfi)

gfi <- gfilinreg(y ~ 0 + group, data = dat, L = 100L, distr = "logistic")
gfiSummary(gfi)

gfi <- gfilinreg(y ~ 0 + group, data = dat, L = 100L, distr = "student", df = 3)
gfiSummary(gfi)


library(heavy)

set.seed(666)
dat <- data.frame(
  group = gl(2, 25),
  y = c(rcauchy(25), 2 + rcauchy(25))
)

(hfit <- heavyLm(y ~ 0 + group, data = dat, family = Cauchy()))
hfit$coefficients
sqrt(hfit$sigma2)

library(MASS)
X <- model.matrix(~ 0 + group, data = dat)
likelihood <- function(y, beta0, beta1, sigma){
  prod(dcauchy((y-X%*%c(beta0,beta1))/sigma)/sigma)
}
(ML <- MASS::fitdistr(dat$y, likelihood, list(beta0=1, beta1=2, sigma=1)))

gfi <- gfilinreg(y ~ 0 + group, data = dat, L = 100L, distr = "cauchy")
gfiSummary(gfi)
