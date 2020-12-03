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
  goodCombs <- matrix(NA_integer_, nrow = p+1L, ncol = 0L)
  for(i in 1:ncol(allCombs)){
    I <- allCombs[, i]
    if(Eigen_rank(X[I, , drop = FALSE]) == p){
      goodCombs <- cbind(goodCombs, I, deparse.level = 0L)
    }
  }
  goodCombs
}

#' @importFrom EigenR Eigen_rank 
#' @noRd
sampleCombinations <- function(X, K){
  n <- nrow(X)
  p <- ncol(X)
  q <- p + 1L
  combs <- matrix(NA_integer_, nrow = p+1L, ncol = 0L)
  combs_chr <- character(0L)
  k <- 0L
  while(k < K){
    I <- sort(sample.int(n, q))
    Ichar <- paste0(I, collapse = "-")
    if(!is.element(Ichar, combs_chr)){
      combs_chr <- c(combs_chr, Ichar)
      if(Eigen_rank(X[I, , drop = FALSE]) == p){
        combs <- cbind(combs, I, deparse.level = 0L)
        k <- k + 1L
      }
    }
  }
  combs
}
