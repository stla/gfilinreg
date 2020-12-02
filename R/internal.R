#' @importFrom Rcpp evalCpp
#' @useDynLib gfilinreg
NULL

#' @importFrom EigenR Eigen_rank 
#' @importFrom utils combn
#' @noRd
goodCombinations <- function(X){
  n <- nrow(X)
  p <- ncol(X)
  allCombs <- combn(n, p+1L)
  goodCombs <- matrix(NA_integer_, nrow = p+1L, ncol = 0)
  for(i in 1:ncol(allCombs)){
    I <- allCombs[, i]
    if(Eigen_rank(X[I, , drop = FALSE]) == p){
      goodCombs <- cbind(goodCombs, I, deparse.level = 0L)
    }
  }
  goodCombs
}
