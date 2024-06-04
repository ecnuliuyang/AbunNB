\name{blackbear}
\alias{blackbear}
\docType{data}
\title{Black bear data}

\description{
This data were collected, over 8 weeks during June and July 2006, to estimate black bear abundance on the Fort Drum Military Installation located in northern New York, USA. Although the experiment was conducted using 38 baited traps, we integrate the encounter histories of 47 individuals from overall traps and treat this data set as a traditional capture--recapture data. Besides the encounter histories, the covariate sex are also available for the bears ever encountered.}

\usage{
  data(blackbear)
}

\format{
  A data frame consisting 47 observations and 9 columns.

\describe{

\item{y1, y2, y3, y4, y5, y6, y7, y8}{
A numeric vector of 0s and 1s for being not captured and captured
on eight capture occasions.}

\item{sex}{
A numeric vector of 0s and 1s for male and female respectively.}
}
}

\details{
This data is used to demonstrate the performance of the penalized empirical likelihood (EL) estimation approach associated with EM algorithms. By fitting a zero-truncated negative binomial regression model to this data, we provide the maximum penalized EL estimator and the EL ratio confidence interval of abundance. The estimation results can be obtained by using the functions \code{\link{abun_nb}}and \code{\link{abun_nb_ci}}.
}

\references{
Gardner, B., Royle, J. A., Wegan, M. T., Rainbolt, R. E., and Curtis, P. D. (2010).
Estimating black bear density using dna data from hair snares. \emph{The Journal of Wildlife
Management}, \bold{74} (2), 318--325.
}

\source{
Royle, J. A., Chandler, R. B., Sollmann, R., and Gardner B. (2014)
\emph{Spatial Capture--Recapture}. Academic Press.
}

\examples{
data(blackbear)
head(blackbear)
table(blackbear$sex)
}

\keyword{datasets}
