library(gfilinreg)

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
