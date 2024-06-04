##' Penalized log EL for the zero-truncated negative binomial regression model
##'
##' @description Function to calculate the penalized empirical log-likelihood (log EL) function for the zero-truncated negative binomial regression model.
##'
##' @param N A number, the population size.
##' @param n A number, the sample size.
##' @param alpha A number, the probability of never being captured.
##' @param prob A vector, the probability masses of covariates.
##' @param numCap A vector, whose elements represent the number of times that individuals are captured.
##' @param z_mat A matrix, where rows represent individuals captured and columns represent covariates.
##' @param beta A vector or matrix with a single column, the regression coefficients in the model \code{"Mh"}.
##' @param k A number, the maximum EL or penalized EL estimates of the dispersion parameter.
##' @param Cp A number, the penalty coefficient when the penalized EL method is used.
##' @param Nchao A number, the Chao (1987)'s lower bound.
##'
##' @return A number, the penalized empirical log-likelihood.
##'
##'
##' @references
##' Chao, A. (1987).
##' Estimating the population size for capture--recapture data with unequal catchability.
##' \emph{Biometrics} \strong{43}, 783-791.
##'
##' @importFrom stats plogis
##'
##' @export
##'
loglikelihood_nb <- function (N, n, alpha, prob, numCap, z_mat, beta, k, Cp, Nchao) {

  y <- as.numeric(numCap)
  zb <- as.numeric( z_mat%*%beta )
  mu <- exp(zb)

  sum( log( N + 1 - c(1:n) ) ) - sum( log(1:n) ) + (N - n)*log( alpha + 1e-300 ) +
    sum( lgamma(y+k) - lgamma(k) - log(factorial(y)) +
                    y * ( log(mu) - log(k+mu) ) + k * ( log(k) - log(k+mu) ) ) +
    sum( log( prob + 1e-300 ) ) -
    Cp * (N - Nchao)^2 * (N > Nchao)

}


##' Maximum penalized EL estimates for the zero-truncated negative binomial regression model
##'
##' @description Function to calculate the maximum penalized empirical likelihood (EL) estimates for the zero-truncated negative binomial regression model using the EM algorithm.
##'
##' @param numCap A vector, whose elements represent the number of times that individuals are captured.
##' @param x A matrix, the individual covariates without the constant one.
##' @param method A character. \code{"EL"} and \code{"PEL"} stand for the EL and penalized EL estimation methods, respectively.
##' @param eps A number, the tolerance of threshold controling the convergence of log-likelihood in EM algorithms.
##' @param maxN A number, the maximum value for \eqn{N} when searching for the estimate of the population size \eqn{N}.
##' @param N0 A number, at which the population size \eqn{N} is fixed when calculating the profile log-likelihood in the \code{\link{abun_nb_ci}} function. It is \code{NULL} when calculating the maximum penalized EL estimates.
##' @param Cp A number, the penalty coefficient when the penalized EL method is used.
##'
##'
##' @return A \code{abun_nb} object.
##'
##'
##' @examples
##'
##' # Fit the zero-truncated negative binomial regression model
##' numCap <- apply(blackbear[,1:8], 1, sum)
##' x <- data.frame(sex = blackbear$sex)
##' bear_el_nb <- abun_nb(numCap = numCap, x = x)
##' bear_pel_nb <- abun_nb(numCap = numCap, x = x, method = "PEL")
##' summary(bear_el_nb)
##' summary(bear_pel_nb)
##'
##' # Calculate the EL ratio confidence interval of the abundance parameter
##' abun_nb_ci(bear_el_nb)
##' abun_nb_ci(bear_pel_nb)
##'
##' @importFrom methods new
##' @importFrom stats coef glm optimize plogis
##' @importFrom MASS glm.nb
##'
##' @export
##'
abun_nb <- function ( numCap, x, method = "EL", eps = 1e-5, maxN = NULL, N0 = NULL, Cp = 0 ) {

  # numCap = object@numCap
  # x = object@x
  # method = object@method
  #
  # eps = object@eps
  # maxN = object@maxN
  # N0 = 300
  # Cp = object@Cp

  ### preparation
  f_fun <- function (y, mu, k) {
    rt <- rep(NA, length(y))
    rt[y == 0] <- k * ( log(k) - log(k+mu[y==0]) )
    rt[y != 0] <- lgamma(y[y!=0]+k) - lgamma(k) - log(factorial(y[y!=0])) +
      y[y!=0] * ( log(mu) - log(k+mu[y!=0]) ) + k * ( log(k) - log(k+mu[y!=0]) )
    rt
  }

  ### initialization
  numCap <- as.numeric(numCap)
  n <- length(numCap)
  if( is.null(maxN) ) maxN <- 100*n


  ### Chao (1987)'s lower bound of N
  f1 <- sum(numCap == 1)
  f2 <- sum(numCap == 2)
  Nchao <- n + f1^2/(2*f2)
  Nchao <- min(Nchao, 1e20)

  if ( method == "PEL" & Cp == 0 )  Cp <- 2 * f2^2 / (n * f1^4)

  x_mat <- as.matrix( x )
  z_mat <- cbind(1, x_mat)
  nx <- rbind(z_mat, z_mat)
  ny <- c(numCap, rep(0,n))

  ### Step 0
  beta <- as.matrix( rep(0, ncol(z_mat)) )
  k <- 1
  prob <- rep(1/n, n)
  zb <- as.numeric(z_mat%*%beta)
  mu <- exp(zb)
  alpha <- sum( exp(f_fun( rep(0,n), mu, k )) * prob )
  N <- ifelse ( is.null(N0), n/( 1 - alpha + 1e-300 ), N0 )

  pars<- c(N, beta, k, alpha)
  likes <- loglikelihood_nb (N, n, alpha, prob, numCap, z_mat, beta, k, Cp, Nchao)

  ### iteration
  err <- 1; nit <- 0

  while ( err > eps ) {

    nit <- nit + 1

    ### calculate omega
    ui <- (N-n) * exp(f_fun( rep(0,n), mu, k )) * prob/( alpha +1e-300 )

    ### update beta
    out <- glm.nb(ny ~ nx - 1, weights = c(rep(1, n), ui), init.theta = k)
    beta <- as.matrix( out$coefficients )
    k <- out$theta
    zb <- as.numeric( z_mat%*%beta )
    mu <- exp(zb)

    ### update prob & alpha
    prob <- (ui + 1) / sum(ui + 1)
    alpha <- sum( exp(f_fun( rep(0,n), mu, k )) * prob )

    ### update N
    if ( is.null(N0) ) {
      obj <- function (nt) sum(log((nt - n + 1):nt)) + (nt - n)*log(alpha + 1e-20) -
        Cp * (nt - Nchao)^2 * (nt > Nchao)
      out <- optimize(obj, lower=n, upper=maxN, maximum=TRUE)
      N <- out$maximum
    }

    ### calculate the log-likelihood
    pars <- rbind(pars, c(N, beta, k, alpha))
    likes <- c(likes, loglikelihood_nb (N, n, alpha, prob, numCap, z_mat, beta, k, Cp, Nchao))

    ### stopping criterion
    err <- likes[nit+1] - likes[nit]

  }

  AIC <- 2*( - likes[nit+1] + 3 + length(beta) )

  rt <- new('abun_nb', N = N, beta = as.numeric(beta), k = k,
            alpha = alpha, method = method,
            loglikelihood = likes[nit+1], AIC = AIC,
            prob = prob, nit = nit, pars = pars, loglikelihoods = likes,
            numCap = numCap, x = x_mat,
            eps = eps, maxN = maxN,
            Cp = Cp, Nchao = Nchao)

  return(rt)

}
