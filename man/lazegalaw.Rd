% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/amen-package.R
\name{lazegalaw}
\alias{lazegalaw}
\title{Lazega's law firm data}
\format{
A list consisting of a socioarray \code{Y} and a nodal attribute matrix \code{X}.

The dyadic variables in \code{Y} include three binary networks: advice, friendship
and co-worker status.

The categorical nodal attributes in \code{X} are coded as follows:
\itemize{
\item status (1=partner, 2=associate)
\item office (1=Boston, 2=Hartford, 3=Providence)
\item practice (1=litigation, 2=corporate)
\item law school (1=Harvard or Yale, 2=UConn, 3=other)
}
\code{seniority} and \code{age} are given in years, and \code{female} is
a binary indicator.
}
\source{
\url{http://moreno.ss.uci.edu/data.html#lazega}
}
\description{
Several nodal and dyadic variables measured on 71 attorneys in a law firm.
}
