#' @name gfilinreg
#' @rdname gfilinreg
#' @title Fiducial sampler for linear regression model
#' @description Weighted samples of the fiducial distribution of the
#'   parameters of a linear regression model with normal, Student, Cauchy, or
#'   logistic error terms.
#'
#' @param formula two-sided formula defining the model
#' @param data dataframe containing the data
#' @param distr the distribution of the error terms, \code{"normal"},
#'   \code{"student"}, \code{"cauchy"}, or \code{"logistic"}
#' @param df degrees of freedom of the Student distribution if
#'   \code{distr = "student"}
#' @param L number of subdivisions of each axis of the hypercube
#'   \code{(0,1)^(p+1)}
#' @param Kmax maximal number of combinations of indices to use
#' @param nthreads number of threads for parallel computations
#' @param stopifbig logical, whether to stop if the algorithm requires huge 
#'   matrices
#'
#' @return A \code{gfilinreg} object, list with the fiducial samples and the
#'   weights.
#'   
#' @references Jan Hannig, Randy C.S. Lai, Thomas C.M. Lee. 
#'   \emph{Computational issues of generalized fiducial inference}.
#'   Computational Statistics and Data Analysis 71 (2014), 849–858.
#'   <doi:10.1016/j.csda.2013.03.003>
#'
#' @examples set.seed(666L)
#' x <- c(1, 2, 3, 4)
#' y <- x + 3 * rcauchy(4L)
#' gfi <- gfilinreg(y ~ x, distr = "cauchy", L = 30L, nthreads = 2L)
#' gfiSummary(gfi)
#'
#' @importFrom EigenR Eigen_rank
#' @importFrom stats model.matrix as.formula
#' @importFrom lazyeval f_eval_lhs f_rhs
#' @importFrom data.table CJ as.data.table rbindlist
#' @importFrom parallel detectCores
#' @export
gfilinreg <- function(
  formula, data = NULL, distr = "student", df = Inf, L = 30L, Kmax = 50L,
  nthreads = parallel::detectCores(), stopifbig = TRUE
){
  if(inSolaris()){
    nthreads <- 1L
  }
  stopifnot(Kmax >= 2)
  distr <- tolower(distr)
  distr <- match.arg(
    distr, c("normal", "gaussian", "student", "cauchy", "logistic")
  )
  if(distr == "student"){
    if(df == Inf){
      distr <- "normal"
    }else if(df == 1){
      distr <- "cauchy"
    }
  }
  y <- f_eval_lhs(formula, data = data)
  X <- model.matrix(formula, data = data)
  betas <- colnames(X)
  X <- unname(X)
  n <- nrow(X)
  p <- ncol(X)
  if(Eigen_rank(X) < p){
    stop("Design is not of full rank.")
  }
  q <- p + 1L
  #
  if(choose(n, q) < 20000L){
    goodCombs <- goodCombinations(X) 
    nGoodCombs <- ncol(goodCombs)
    if(Kmax >= nGoodCombs){
      K <- nGoodCombs
      combs <- goodCombs
    }else{
      combs <- goodCombs[, sample.int(nGoodCombs, Kmax)]
      K <- Kmax
    }
  }else{
    message(sprintf("Sampling %d combinations of indices...", Kmax))
    combs <- sampleCombinations(X, Kmax)
    message("Done.")
    K <- Kmax
  }
  #
  if(!enoughRAM(q*L^q*(1+K/2))){
    stop("Insufficient free RAM.")
  }
  #
  if(stopifbig && (K * q * L^q / 2 > 8e7)){
    stop(
      paste0(
        "The algorithm needs to deal with two big matrices ",
        sprintf(
          "(the biggest with %g entries: %g GB). ", 
          K * q * L^q / 2, K * q * L^q / 2 * 8e-9
        ), 
        "Set the option `stopifbig` to FALSE if you want to proceed anyway."
      )
    )
  }
  M <- ceiling(L^q / 2L)
  # centers of hypercubes (volume 1/L^p)
  centers <- as.matrix(
    do.call(
      function(...){CJ(..., sorted = FALSE)}, 
      rep(list(seq(0, 1, length.out = L+1L)[-1L] - 1/(2*L)), q)
    )
  )
  # algorithm
  XIs <- do.call(cbind, lapply(1L:K, function(k) X[combs[, k], , drop = FALSE]))
  XmIs <- 
    do.call(cbind, lapply(1L:K, function(k) X[-combs[, k], , drop = FALSE]))
  yIs <- apply(combs, 2L, function(I) y[I])
  ymIs <- rbind(apply(combs, 2L, function(I) y[-I]))
  if(distr %in% c("normal", "gaussian")){
    cpp <- f_normal(
      centers = t(centers),
      XIs = XIs, XmIs = XmIs,
      yIs = yIs, ymIs = ymIs,
      K = K, p = p, M = M, n = n, 
      nthreads = nthreads
    )
  }else if(distr == "student"){
    cpp <- f_student(
      centers = t(centers),
      XIs = XIs, XmIs = XmIs,
      yIs = yIs, ymIs = ymIs,
      K = K, p = p, M = M, n = n,
      nu = df,
      nthreads = nthreads
    )
  }else if(distr == "cauchy"){
    cpp <- f_cauchy(
      centers = t(centers),
      XIs = XIs, XmIs = XmIs,
      yIs = yIs, ymIs = ymIs,
      K = K, p = p, M = M, n = n,
      nthreads = nthreads
    )
  }else if(distr == "logistic"){
    cpp <- f_logistic(
      centers = t(centers),
      XIs = XIs, XmIs = XmIs,
      yIs = yIs, ymIs = ymIs,
      K = K, p = p, M = M, n = n,
      nthreads = nthreads
    )
  }
  #
  LOGWEIGHTS <- lapply(cpp, `[[`, "logWeights")
  THETAS <- lapply(cpp, function(output){
    as.data.table(`colnames<-`(output[["Theta"]], c(betas, "sigma")))
  })
  J <- exp(do.call(c, LOGWEIGHTS))
  out <- list(
    Theta = rbindlist(THETAS),
    weight = J/sum(J)
  )
  attr(out, "distr") <- distr
  attr(out, "df") <- df
  rhs <- as.character(f_rhs(formula))
  if(rhs[1L] == "+") rhs <- rhs[-1L]
  attr(out, "formula") <- as.formula(
    paste0("~ ", paste0(rhs, collapse = " + "))
  )
  class(out) <- "gfilinreg"
  out
}
