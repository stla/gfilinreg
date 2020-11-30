inference <- function(fidsamples, param, alpha = 0.05){
  fipred <- inherits(fidsamples, "gfilinreg.pred")
  out <- numeric(4L)
  names(out) <- c("mean", "median", "lwr", "upr")
  sample <-
    if(fipred) fidsamples[["FPD"]][[param]] else fidsamples[["Theta"]][[param]]
  weights <- fidsamples[["weight"]]
  out[1L] <- sum(sample * weights) # mean
  h <- cbind(sample, weights)
  hsort <- h[order(h[,1L]), ]
  hsum <- cumsum(hsort[, 2L])
  ci_u <- min(which(hsum >= 1-alpha/2))
  ci_l <- min(which(hsum >= alpha/2))
  ci_m <- min(which(hsum >= 0.5))
  out[3L] <- hsort[ci_l, 1L] # lower bound
  out[4L] <- hsort[ci_u, 1L] # upper bound
  out[2L] <- hsort[ci_m, 1L] # estimate (median)
  out
}

#' Title
#'
#' @param fidsamples xx
#' @param conf confidence level
#'
#' @return xx
#' @export
#'
#' @examples xx
gfiSummary <- function(fidsamples, conf = 0.95){
  sims <- if(inherits(fidsamples, "gfilinreg.pred")){
    fidsamples[["FPD"]]
  }else fidsamples[["Theta"]]
  seq_ <- 1L:ncol(sims)
  names(seq_) <- names(sims)
  out <-
    t(vapply(seq_, function(x) inference(fidsamples, x, 1-conf), numeric(4L)))
  attr(out, "confidence level") <- conf
  out
}

#' Title
#'
#' @param parameter xx
#' @param fidsamples xx
#'
#' @return xx
#'
#' @importFrom lazyeval f_eval_rhs
#' @importFrom spatstat ewcdf
#' @export
#'
#' @examples xx
gfiCDF <- function(parameter, fidsamples){
  dataName <- ifelse(inherits(fidsamples, "gfilinreg.pred"), "FPD", "Theta")
  data <- fidsamples[[dataName]]
  fsims <- f_eval_rhs(parameter, data = data)
  ewcdf(fsims, weights = fidsamples[["weight"]])
}

#' Title
#'
#' @param parameter xx
#' @param fidsamples xx
#' @param conf xx
#'
#' @return xx
#'
#' @importFrom spatstat quantile.ewcdf
#' @export
#'
#' @examples xx
gfiConfInt <- function(parameter, fidsamples, conf = 0.95){
  fcdf <- fiCDF(parameter, fidsamples)
  alpha <- 1 - conf
  quantile.ewcdf(fcdf, c(alpha/2, 1-alpha/2))
}

#' Title
#'
#' @param parameter xx
#' @param fidsamples xx
#' @param probs xx
#'
#' @return xx
#'
#' @importFrom spatstat quantile.ewcdf
#' @export
#'
#' @examples xx
gfiQuantile <- function(parameter, fidsamples, probs){
  fcdf <- gfiCDF(parameter, fidsamples)
  quantile.ewcdf(fcdf, probs = probs)
}

