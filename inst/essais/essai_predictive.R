library(filinreg)

set.seed(666)
x <- rnorm(15)
y <- x + rnorm(15)
new <- data.frame(x = c(-3,3))
predict(lm(y ~ x), new, interval = "prediction")


fi <- filinreg(y ~ x, L = 100L, distr = "normal", lucky = TRUE)

fpred <- filinregPredictive(fi, new)
fiSummary(fpred)

