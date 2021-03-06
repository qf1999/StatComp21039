% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gibbsC}
\alias{gibbsC}
\title{A Gibbs sampler using Rcpp}
\usage{
gibbsC(N, thin, a, b, n)
}
\arguments{
\item{N}{the number of samples}

\item{thin}{the number of between-sample random numbers}

\item{a}{the parameter for gibbsc}

\item{b}{the parameter for gibbsc}

\item{n}{the parameter for gibbsc}
}
\value{
a random sample from gibbs \code{mat}
}
\description{
A Gibbs sampler using Rcpp.
}
\examples{
\dontrun{
    gc <- gibbsC(100, 10, 2, 4, 16)
    print(gc)
}
}
