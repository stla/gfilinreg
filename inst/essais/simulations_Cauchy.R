library(gfilinreg)
library(heavy)
library(data.table)

set.seed(666L)
nsims <- 500L

MLD <- matrix(NA_real_, nrow = nsims, ncol = 3L)
colnames(MLD) <- c("group1", "group2", "sigma")
FIDlist <- vector("list", length = nsims)

for(i in 1L:nsims){
  cat(i, " - ")
  dat <- data.frame(
    group = gl(2, 6),
    y = c(rcauchy(6), 2 + rcauchy(6))
  )
  #
  hfit <- heavyLm(y ~ 0 + group, data = dat, family = Cauchy())
  MLD[i, ] <- c(hfit$coefficients, sqrt(hfit$sigma2))
  #
  gfi <- gfilinreg(y ~ 0 + group, data = dat, L = 100L, distr = "cauchy")
  FIDlist[[i]] <- cbind(
    param = c("group1", "group2", "sigma"), 
    as.data.table(gfiSummary(gfi))
  )
}

FID <- rbindlist(FIDlist)

saveRDS(MLD, "~/Work/R/gfilinreg/inst/essais/MLD.rds")
saveRDS(FID, "~/Work/R/gfilinreg/inst/essais/FID.rds")
saveRDS(FIDlist, "~/Work/R/gfilinreg/inst/essais/FIDlist.rds")

# FID[param=="group1"] etc

