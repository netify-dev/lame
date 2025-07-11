% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/amen-package.R
\docType{data}
\name{YX_bin}
\alias{YX_bin}
\title{binary relational data and covariates}
\format{
The format is: List of 2 $ Y: num \link{1:100, 1:100} NA 0 0 0 0 0 0 0 0
1 ...  $ X: num \link{1:100, 1:100, 1:8} 1 1 1 1 1 1 1 1 1 1 ...  ..- attr(*,
"dimnames")=List of 3 .. ..$ : NULL .. ..$ : NULL .. ..$ : chr \link{1:8}
"intercept" "rgpa" "rsmoke" "cgpa" ...
}
\usage{
data(YX_bin)
}
\description{
a synthetic dataset that includes binary relational data as well as
information on eight covariates
}
\examples{

data(YX_bin)
gofstats(YX_bin$Y) 

}
\keyword{datasets}
