% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BiLoglik.R
\name{LL.BiInvPa}
\alias{LL.BiInvPa}
\title{Log-likelihood function for bivariate inverted-Pareto density}
\usage{
LL.BiInvPa(par, x, family)
}
\arguments{
\item{par}{Numeric; a parameter for bivariate inverted-Pareto densities}

\item{x}{A \eqn{N\times 2} data matrix from a bivariate inverted-Pareto distribution}

\item{family}{An integer; indicates the type of bivariate tail copula families.
Possible tail copula families include:
1=Husler-Reiss
2=Negative Logistic
3=Logistic
4=Dirichlet}
}
\value{
Numeric; a negative log-likelihood value
}
\description{
\code{LL.BiInvPa} calculates a negative log-likelihood of a bivariate inverted-Pareto density.
For details on the inverted-Pareto distribution, refer to Kiriliouk, A., Lee, J., & Segers, J. (2023).
}
\references{
Kiriliouk, A., Lee, J., & Segers, J. (2023). X-Vine Models for Multivariate Extremes. arXiv preprint arXiv:2312.15205.
}
