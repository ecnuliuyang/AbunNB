##' abun_nb class
##'
##' @description A S4 class is used to present the data and maximum penalized empirical likelihood (EL) estimates under the zero-truncated negative binomial regression model.
##' @slot method A character. \code{"EL"} and \code{"PEL"} stand for the EL and penalized EL estimation methods, respectively.
##' @slot N A number, the maximum EL or PEL estimate of the population size.
##' @slot Nchao A number, the Chao (1987)'s lower bound.
##' @slot beta A vector, the maximum EL or penalized EL estimates of regression coefficients.
##' @slot k A number, the maximum EL or penalized EL estimates of the dispersion parameter.
##' @slot alpha A number, the maximum EL or penalized EL estimate of the probability of never being captured.
##' @slot loglikelihood A number, the log EL or penalized log EL value.
##' @slot AIC A number, the Akaike Information Criterion value.
##' @slot prob A vector, the probability masses of individual covariates.
##' @slot nit A number, the number of iterations of EM algorithms.
##' @slot pars A matrix. Row represent the iterative values of parameters.
##' @slot loglikelihoods A vector, the log-likelihood values in iteration.
##' @slot numCap A vector, whose elements represent the number of times that individuals are captured.
##' @slot x A matrix, the individual covariates without the constant one.
##' @slot eps A number, the tolerance of threshold controling the convergence of log-likelihood in EM algorithms.
##' @slot maxN A number, the maximum value for \eqn{N} when searching for the estimate of the population size \eqn{N}.
##' @slot Cp A number, the penalty coefficient when the penalized EL method is used.
##'
##' @export
##'
##' @references
##' Chao, A. (1987).
##' Estimating the population size for capture--recapture data with unequal catchability.
##' \emph{Biometrics} \strong{43}, 783-791.
##'
setClass( 'abun_nb', slots =
            list( method = "character",
                  N = "numeric",
                  Nchao = "numeric",
                  beta = "numeric",
                  k = "numeric",
                  alpha = "numeric",
                  loglikelihood = "numeric",
                  AIC =  "numeric",
                  prob = "numeric",
                  nit = "numeric",
                  pars = "matrix",
                  loglikelihoods = "numeric",
                  numCap = "numeric",
                  x = "matrix",
                  eps = "numeric",
                  maxN = "numeric",
                  Cp = "numeric" ) )

##' Function to show the \code{abun_nb} object
##'
##' @param object \code{abun_nb} object
setMethod("show", "abun_nb",

          function(object){

            coefficients <- round(object@beta, 4L)
            cat(paste0("\nMaximum ", toupper(object@method)),
                "estimation for the zero-truncated negative binomial regression model\n\n")

            cat ("Maximum", object@method, "estimate of N:", round(object@N))
            cat ("\n")
            cat ("Maximum", object@method, "estimate of beta:", coefficients)
            cat ("\n")
            cat ("Maximum", object@method, "estimate of k:", round(object@k, 4L))
            cat ("\n\n")
            cat ("Log-likelihood:", round(object@loglikelihood, 2L))
            cat ("\n")
            cat ("AIC:", round(object@AIC, 2L))

          })


##' abun_nb_summary class
##'
##' @description A S4 class summarizes the maximum (penalized) empirical likelihood estimates for the zero-truncated negative binomial regression model
##' @slot coefficients A matrix, which summarizes the estimates of regression coefficients.
##' @slot dispersion A matrix, which summarizes the estimates of the dispersion parameter.
##' @slot abundance A matrix, which summarizes the estimates of the population size.
##'
##' @export
##'
setClass( "abun_nb_summary", contains = "abun_nb",
          slots = list( coefficients = "matrix",
                        dispersion = "matrix",
                        abundance = "matrix" ) )


##' Function to summarize the maximum (penalized) empirical likelihood estimates for the zero-truncated negative binomial regression model
##'
##' @param object a \code{abun_nb} object
##'
##' @importFrom methods new
##' @importFrom stats pnorm
##'
##' @export
##'
setMethod('summary', 'abun_nb',
          function(object) {
            rt <- new('abun_nb_summary',
                      method = object@method,
                      N = object@N,
                      Nchao = object@Nchao,
                      beta = object@beta,
                      k = object@k,
                      alpha = object@alpha,
                      loglikelihood = object@loglikelihood,
                      AIC = object@AIC,
                      prob = object@prob,
                      nit = object@nit,
                      pars = object@pars,
                      loglikelihoods = object@loglikelihoods,
                      numCap = object@numCap,
                      x = object@x,
                      eps = object@eps,
                      maxN = object@maxN,
                      Cp = object@Cp)

            est_beta <- object@beta
            se <- abun_nb_se(object)
            se_beta <- se$se_beta
            zval <- est_beta/abs(se_beta)
            coefficients <- cbind(est_beta, se_beta, zval, 2L * pnorm(abs(zval), lower.tail = FALSE))
            colnames(coefficients) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')

            name_x <- names(as.data.frame(rt@x))
            rownames(coefficients) <- c( "Intercept", name_x )

            rt@coefficients <- coefficients

            dispersion <- matrix(c(round(object@k,2), round(se$se_k,2L)), 1L)
            colnames(dispersion) <- c("Estimate", "Std. Error" )
            rownames(dispersion) <- '      '
            rt@dispersion <- dispersion

            abundance <- matrix(c(round(object@N), round(se$se_N,2L)), 1L)
            colnames(abundance) <- c("Estimate", "Std. Error" )
            rownames(abundance) <- '      '
            rt@abundance <- abundance
            return(rt)
          })


##' Function to show the \code{abun_nb_summary} object
##'
##' @param object \code{abun_nb_summary} object
##' @importFrom stats printCoefmat
##'
setMethod('show', 'abun_nb_summary',
          function(object) {
            digits = max(3L, getOption("digits") - 3L)
            signif.stars = getOption("show.signif.stars")

            cat(paste0("\nMaximum ", toupper(object@method)),
                "estimates for the zero-truncated negative binomial regression model\n")

            cat("Coefficients:\n")
            printCoefmat(object@coefficients, digits = digits, signif.stars = signif.stars, na.print = "NA")

            cat("\nDispersion:\n")
            print(object@dispersion)

            cat("\nAbundance:\n")
            print(object@abundance)

            cat ('\nLog-likelihood:', object@loglikelihood,
                 'with the AIC value being', round(object@AIC, 2L),
                 '\n\n')
            invisible(object)
          })
