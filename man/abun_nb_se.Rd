% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abun_nb_se.R
\name{abun_nb_se}
\alias{abun_nb_se}
\title{Standard errors for the zero-truncated negative binomial regression model}
\usage{
abun_nb_se(object)
}
\arguments{
\item{object}{A \code{abun_nb} object.}
}
\value{
A list with three elements:

\itemize{
  \item se_N, the standard error for the population size.
  \item se_beta, the standard error for the regression coefficients.
  \item se_k, the standard error for the dispersion parameter.
}
}
\description{
Function to calculate the standard errors, based on the asymptotic normality, of the maximum penalized empirical likelihood estimates for the zero-truncated negative binomial regression model
}
