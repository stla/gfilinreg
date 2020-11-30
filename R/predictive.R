#' Title
#'
#' @param fidsamples xx
#' @param newdata xx
#'
#' @return xx
#'
#' @importFrom stats model.matrix rt rlogis rcauchy rnorm
#' @export
#'
#' @examples xx
gfilinregPredictive <- function(fidsamples, newdata){
  if(is.null(newdata) || missing(newdata)){
    newdata <- as.data.frame(matrix(nrow = 1L, ncol = 0L))
  }
  if(anyDuplicated(newdata)){
    stop(
      "There are some duplicated rows in `newdata`."
    )
  }
  # newdata <- droplevels(newdata)
  X <- model.matrix(attr(fidsamples, "formula"), data = newdata)
  N <- length(fidsamples[["weight"]])
  Theta <- t(as.matrix(fidsamples[["Theta"]]))
  Beta <- Theta[-nrow(Theta), ]
  sigma <- Theta[nrow(Theta), ]
  if(attr(fidsamples, "distr") == "student"){
    rdistr <- function(n) rt(n, df = attr(fidsamples, "df"))
  }else if(attr(fidsamples, "distr") == "logistic"){
    rdistr <- function(n) rlogis(n)
  }else if(attr(fidsamples, "distr") == "cauchy"){
    rdistr <- function(n) rcauchy(n)
  }else{
    rdistr <- function(n) rnorm(n)
  }
  n <- nrow(newdata)
  out <- matrix(NA_real_, nrow = N, ncol = n)
  colnames(out) <- paste0("y", seq_len(n))
  for(i in 1L:N){
    out[i, ] <- X %*% Beta[, i] + sigma[i]*rdistr(n)
  }
  out <- list(FPD = as.data.frame(out), weight = fidsamples[["weight"]])
  class(out) <- c("gfilinreg", "gfilinreg.pred")
  out
}
