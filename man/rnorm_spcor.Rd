% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_functions.R
\name{rnorm_spcor}
\alias{rnorm_spcor}
\title{Generate spatially-autocorrelated random variables}
\usage{
rnorm_spcor(n, covar, fast = FALSE)
}
\arguments{
\item{n}{number of samples per location}

\item{covar}{a matrix containing pairwise covariances among
locations.}

\item{fast}{logical: should the faster method, whose realizations
observe lower fidelity to the covariance function than when
\code{fast = FALSE}.}
}
\value{
a matrix with rows corresponding to locations and columns
corresponding to random samples taken at that location
}
\description{
Generate spatially-autocorrelated random variables
}
\examples{
cor.mat <- matrix(c(1.0, 0.3, 0.1,
                    0.3, 1.0, 0.0,
                    0.1, 0.0, 1.0), nrow = 3)
rnorm_spcor(5, cor.mat)

## plot sample covariance vs. true covariance
spcov = generate_spcov(c(10, 10), sd = 1, unit.map = TRUE)
D = spcov$D.mat
C = spcov$covar.mat
x <- rnorm_spcor(100, C)
## sample covariance vs. distance
plot(x = D, y = cor(t(x)), xlab = "distance", ylab = "covariance")
## add line for covariance function
curve(remotePARTS::covar_exp(x, 1), from = 0, to = max(D),
      add = TRUE, col = "red")

}
