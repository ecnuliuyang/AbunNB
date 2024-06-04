##' Standard errors for the zero-truncated negative binomial regression model
##'
##' @description Function to calculate the standard errors, based on the asymptotic normality, of the maximum penalized empirical likelihood estimates for the zero-truncated negative binomial regression model
##'
##' @param object A \code{abun_nb} object.
##'
##' @return A list with three elements:
##'
##' \itemize{
##'   \item se_N, the standard error for the population size.
##'   \item se_beta, the standard error for the regression coefficients.
##'   \item se_k, the standard error for the dispersion parameter.
##' }
##'
##' @export
##'

abun_nb_se <- function (object) {

  ### preparation
  f_fun <- function (y, xb) {
    mu <- exp(xb)
    unlist(lapply(y, function(x) factorial(k + x - 1)/factorial(x)/factorial(k - 1))) * (k/(k + mu))^k * (mu/(k + mu))^y
  }


  s2_fun <- function(y) {
    if (y == 0) {
      out <- 0
    } else {
      c <- numeric(y)
      for (i in 0:(y-1)) {
        c[i + 1] <- (k + i)^(-2)
      }
      out <- -sum(c)
    }
    out
  }

  y = object@numCap
  x_mat <- object@x
  x_mat <- as.matrix( x_mat[order(y),] )
  z_mat <- cbind( 1, x_mat )
  y <- sort(y)
  n <- length(y)
  N <- object@N
  beta <- as.matrix(object@beta)
  k <- object@k
  alp <- object@alpha

  zb <- as.numeric(z_mat%*%beta)
  s <- length(beta)

  f0 <- f_fun(rep(0,n), zb)
  mu <- exp(zb)

  varphi <- sum((1-f0)^(-2))/N
  psi <- log(k/(k+mu)) + mu/(k+mu)
  s2 <- unlist(lapply(y,s2_fun))
  zr <- rep(0, s)
  zr <- matrix(zr, nrow = s, ncol = 1)

  ### Vij's
  V11 <- 1-1/alp
  V41 <- 1/alp
  V22 <- t(z_mat) %*% diag(
    (k*f0*mu/(1-f0)-k-mu)*k*mu/(k+mu)^2/(1-f0) , n, n ) %*% z_mat/N

  V23 <- -t(z_mat)%*%(k*f0*psi*mu/(k+mu)/(1-f0)^2)/N
  V24 <- t(z_mat)%*%(f0*k*mu/(k+mu)/(1-f0)^2)/N
  V25 <- (1-alp)^2*V24

  V33 <- sum((f0*psi^2/(1-f0) + (k*mu+mu^2)/k/(k+mu)^2)/(1-f0))/N + sum(s2/(1-f0))/N
  V34 <- - sum(f0*psi/(1-f0)^2)/N
  V35 <- (1-alp)^2*V34

  V44 <- varphi -1/alp
  V45 <- (1-alp)^2*varphi
  V55 <- (1-alp)^4*varphi - (1-alp)^3

  V55i <- 1/V55

  W <- rbind(cbind(-V11,t(zr),0,-V41),
             cbind(zr, -V22+V55i*V25%*%t(V25), -V23+V55i*V35*V25, -V24+V55i*V45*V25),
             cbind(0, -t(V23)+V35*V55i*t(V25),-V33+V35*V55i*V35,-V34+V35*V55i*V45),
             cbind(-V41,-t(V24)+V45*V55i*t(V25),-V34+V45*V55i*V35,-V44+V45*V55i*V45))

  se <- diag(solve(W))
  se_N <- sqrt(se[1]*N)
  se_beta <- sqrt(se[2:(s+1)]/N)
  se_k <- sqrt(se[(s+2)]/N)
  rt <- list( se_N = se_N, se_beta = se_beta, se_k=se_k)
  return(rt)
}
