library(gfilinreg)

set.seed(666)
x <- rnorm(15)
y <- x + rnorm(15)

gfi <- gfilinreg(y ~ x, L = 20L, distr = "normal", lucky = TRUE)
gfiSummary(gfi)
gfiR <- gfilinreg:::gfilinregR(y ~ x, L = 20L, distr = "student", df = Inf, lucky = TRUE)
gfiSummary(gfiR) # idem


gfi <- gfilinreg(y ~ x, L = 20L, distr = "student", df = 3, lucky = TRUE)
gfiSummary(gfi)
gfiR <- gfilinreg:::gfilinregR(y ~ x, L = 20L, distr = "student", df = 3, lucky = TRUE)
gfiSummary(gfiR) # idem


gfi <- gfilinreg(y ~ x, L = 20L, distr = "cauchy", lucky = TRUE)
gfiSummary(gfi)
gfiR <- gfilinreg:::gfilinregR(y ~ x, L = 20L, distr = "student", df = 1, lucky = TRUE)
gfiSummary(gfiR) # idem


gfi <- gfilinreg(y ~ x, L = 20L, distr = "logistic", lucky = TRUE)
gfiSummary(gfi)
gfiR <- gfilinreg:::gfilinregR(y ~ x, L = 20L, distr = "logistic", lucky = TRUE)
gfiSummary(gfiR) # idem
