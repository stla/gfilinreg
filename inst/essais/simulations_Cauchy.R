library(gfilinreg)
library(heavy)
library(data.table)

nsims <- 500L
MAXLHD <- matrix(NA_real_, nrow = nsims, ncol = 3L)
colnames(MAXLHD) <- c("group1", "group2", "sigma")
FIDlist <- vector("list", length = nsims)

group <- gl(2L, 3L)
set.seed(666L)
for(i in 1L:48){
  cat(i, " - ")
  # simulated dataset
  dat <- data.frame(
    group = group,
    y = c(rcauchy(3L), 2 + rcauchy(3L))
  )
  # # max-likelihood estimates
  # hfit <- heavyLm(y ~ 0 + group, data = dat, family = Cauchy())
  # MAXLHD[i, ] <- c(hfit[["coefficients"]], sqrt(hfit[["sigma2"]]))
  # # fiducial stuff
  # fidsamples <- gfilinreg(y ~ 0 + group, data = dat, L = 50L, distr = "cauchy")
  # FIDlist[[i]] <- cbind(
  #   parameter = c("group1", "group2", "sigma"), 
  #   as.data.table(gfiSummary(fidsamples))
  # )
  # rm(list = "fidsamples")
}

fidsamples <- gfilinreg(y ~ 0 + group, data = dat, L = 50L, distr = "cauchy")

stop() 


FID <- rbindlist(FIDlist)

saveRDS(MAXLHD, "~/Work/R/gfilinreg/inst/essais/MAXLHD.rds")
saveRDS(FID, "~/Work/R/gfilinreg/inst/essais/FID.rds")
saveRDS(FIDlist, "~/Work/R/gfilinreg/inst/essais/FIDlist.rds")

stop()


MAXLHD <- readRDS("~/Work/R/gfilinreg/inst/essais/MAXLHD.rds")
FID <- readRDS("~/Work/R/gfilinreg/inst/essais/FID.rds")

#usethis::use_data(MAXLHD, FID)

################################################################################
library(data.table)
data("FID")
data("MAXLHD")

library(kde1d)
group1_maxlhd     <- MAXLHD[, "group1"]
group1_fid_mean   <- FID[parameter == "group1"][["mean"]]
group1_fid_median <- FID[parameter == "group1"][["median"]]

kfit_maxlhd     <- kde1d(group1_maxlhd, mult = 4)
kfit_fid_mean   <- kde1d(group1_fid_mean, mult = 4)
kfit_fid_median <- kde1d(group1_fid_median, mult = 4)

curve(
  dkde1d(x, kfit_maxlhd), from = -4, to = 4, 
  lwd = 2, col = "red", lty = "dashed"
)
curve(
  dkde1d(x, kfit_fid_mean), add = TRUE, 
  lwd = 2, col = "green", lty = "dashed"
)
curve(
  dkde1d(x, kfit_fid_median), add = TRUE, 
  lwd = 2, col = "blue", lty = "dashed"
)


