% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{thomas_algorithm}
\alias{thomas_algorithm}
\title{Thomas algorithm for solving tridiagonal linear systems. Stability only for diagonally dominant or symmetric positive definite.}
\usage{
thomas_algorithm(a, b, c, d)
}
\arguments{
\item{a}{lower diagonal of matrix}

\item{b}{diagonal of matrix}

\item{c}{upper diagonal of matrix}

\item{d}{vector on RHS of equation}
}
\description{
{Thomas algorithm for tridiagonal linear systems. This can be used to solve
PDEs arising from Feynman-Kac connections.}
}
\details{
{The algorithm is on wikipedia and quite straightforward.}
}
