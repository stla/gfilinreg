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
#' @param lucky logical, whether to perform the matrix inversions in the
#'   algorithm without checking invertibility; it is possible that some of
#'   these matrices are not invertible, if you are unlucky
#'
#' @return A \code{gfilinreg} object, list with the fiducial samples and the
#'   weights.
#'
#' @examples set.seed(666)
#' x <- c(1, 2, 3, 4)
#' y <- x + 3 * rcauchy(4L)
#' gfi <- gfilinreg(y ~ x, distr = "cauchy", L = 30L)
#' gfiSummary(gfi)
#'
#' @importFrom arrangements icombinations
#' @importFrom EigenR Eigen_rank
#' @importFrom stats model.matrix as.formula
#' @importFrom lazyeval f_eval_lhs f_rhs
#' @importFrom data.table CJ as.data.table
#' @export
gfilinreg <- function(
  formula, data = NULL, distr = "student", df = Inf, L = 10L, lucky = TRUE
){
  distr <- match.arg(distr, c("normal", "student", "cauchy", "logistic"))
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
  # centers of hypercubes (volume 1/L^p)
  centers <- as.matrix(
    do.call(
      CJ, rep(list(seq(0, 1, length.out = L+1L)[-1L] - 1/(2*L)), q)
    )
  )
  # select indices
  Iiterator <- icombinations(n, q)
  I <- Iiterator$getnext()
  XI <- X[I, , drop = FALSE]
  while(Eigen_rank(XI) < p){
    I <- Iiterator$getnext()
    XI <- X[I, , drop = FALSE]
  }
  XmI <- X[-I, , drop = FALSE]
  yI <- y[I]
  ymI <- y[-I]
  # algorithm
  if(Eigen_rank(cbind(XI, 1)) < q){
    # remove centers having equal coordinates (H'H is not invertible)
    centers <-
      centers[apply(centers, 1L, function(row) length(unique(row)) > 1L),]
    M <- (L^q - L) / 2L # number of centers yielding sigma>0
  }else{
    M <- floor(L^q / 2L) # TODO: test !!! - done, seems OK
  }
  if(lucky){
    if(distr == "normal"){
      cpp <- f_normal(
        centers = t(centers),
        XI = XI, XmI = XmI,
        yI = yI, ymI = ymI,
        M = M, n = n
      )
    }else if(distr == "student"){
      cpp <- f_student(
        centers = t(centers),
        XI = XI, XmI = XmI,
        yI = yI, ymI = ymI,
        M = M, n = n,
        nu = df
      )
    }else if(distr == "cauchy"){
      cpp <- f_cauchy(
        centers = t(centers),
        XI = XI, XmI = XmI,
        yI = yI, ymI = ymI,
        M = M, n = n
      )
    }else if(distr == "logistic"){
      cpp <- f_logistic(
        centers = t(centers),
        XI = XI, XmI = XmI,
        yI = yI, ymI = ymI,
        M = M, n = n
      )
    }
  }else{
    # xxx
  }
  J <- exp(cpp[["logWeights"]])
  out <- list(
    Theta = as.data.table(`colnames<-`(cpp[["Theta"]], c(betas, "sigma"))),
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
