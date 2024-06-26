% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ParRange.R
\name{ParRangeTC}
\alias{ParRangeTC}
\title{Range of parameters for tail copula families}
\usage{
ParRangeTC(family)
}
\arguments{
\item{family}{Integer; indicates the type of bivariate tail copula densities.
Possible tail copula families include:
1=Husler-Reiss
2=Negative Logistic
3=Logistic
4=Dirichlet}
}
\value{
A numeric vector; specifies the range of parameter for the corresponding tail copula density
}
\description{
\code{ParRangeTC} specifies a proper range of parameters for tail copula families.
This range specification helps the convergence of the \code{Optimize} command.
}
