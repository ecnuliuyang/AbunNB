##' Penalized EL ratio confidence interval for the zero-truncated negative binomial regression model
##'
##' @description Penalized empirical likelihood (EL) ratio confidence interval of the abundance parameter \eqn{N} is calculated for the zero-truncated negative binomial regression model.
##'
##' @param object A \code{abun_nb} object.
##' @param level A number, the nominal level. Default is 0.95.
##' @param UB A number, the upper bound when searching for the upper limit of \eqn{N}. Default is \eqn{1e+09}.
##'
##' @return A vector, the penalized EL ratio confidence interval of \eqn{N}.
##'
##' @importFrom stats qnorm qchisq uniroot
##' @importFrom MASS glm.nb
##'
##' @export
##'
abun_nb_ci <- function ( object, level = 0.95, UB = 1e9 ) {

  ###  maximum penalized empirical log-likelihood
  like_full <- object@loglikelihood
  n <- length(object@numCap)

  rn <- function (N0) {

    like_null <- abun_nb(numCap = object@numCap, x = object@x, method = object@method,
                         eps = object@eps, maxN = object@maxN, N0 = N0,
                         Cp = object@Cp)

    2 * ( like_full - like_null@loglikelihood ) - qchisq(level, 1)

  }

  hatN <- object@N
  ntemp <- c(hatN, 10 * hatN)
  nit <- 2
  ind <- rn(ntemp[nit]) <= 0
  while ( ind ) {

    ntemp <- c(ntemp, ntemp[nit]*10)
    nit <- nit+1
    ind <- rn(ntemp[nit]) <= 0

    if (ntemp[nit] >= UB) {
      ci_upper = UB
      break
    }

  }

  if (!ind & ntemp[nit] < UB)  ci_upper <- uniroot( rn, c(ntemp[nit-1], ntemp[nit]), tol=0.01 )$root

  if ( rn(n) <= 0 ) {
    ci_lower <- n
  } else {
    ci_lower <- uniroot( rn, c(n, hatN), tol=0.01 )$root
  }

  return( c(ci_lower, ci_upper) )

}



