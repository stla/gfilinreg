library(gfilinreg)

set.seed(666)
x <- rnorm(15)
y <- x + rnorm(15)
new <- data.frame(x = c(-3,3))
predict(lm(y ~ x), new, interval = "prediction")


gfi <- gfilinreg(y ~ x, L = 100L, distr = "normal")

fpred <- gfilinregPredictive(gfi, new)
gfiSummary(fpred)

