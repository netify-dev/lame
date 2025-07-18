% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/amen-package.R
\docType{data}
\name{YX_nrm}
\alias{YX_nrm}
\title{normal relational data and covariates}
\format{
The format is: List of 2 $ Y: num \link{1:100, 1:100} NA -4.05 -0.181
-3.053 -1.579 ...  $ X: num \link{1:100, 1:100, 1:8} 1 1 1 1 1 1 1 1 1 1 ...  ..-
attr(*, "dimnames")=List of 3 .. ..$ : NULL .. ..$ : NULL .. ..$ : chr \link{1:8}
"intercept" "rgpa" "rsmoke" "cgpa" ...
}
\usage{
data(YX_nrm)
}
\description{
a synthetic dataset that includes continuous (normal) relational data as
well as information on eight covariates
}
\examples{

data(YX_nrm)
gofstats(YX_nrm$Y)


}
\keyword{datasets}
