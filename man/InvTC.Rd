% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BiTCfamily.R
\name{InvTC}
\alias{InvTC}
\title{Quantile function of a bivariate tail copula}
\usage{
InvTC(w2, u1, par, family)
}
\arguments{
\item{w2}{A numeric vector of the first argument for the bivariate tail copula.}

\item{u1}{A numeric vector of the second argument conditioned on \eqn{u_1 < 1}
for the bivariate tail copula.}

\item{par}{A numeric of parameter values.}

\item{family}{A numeric vector indicating the list of tail copula families.}
}
\value{
A numeric vector of quantile values.
}
\description{
Conditioning on the second argument of the bivariate tail copula being less than 1,
\eqn{x_1 < 1}, the function \code{InvTC()} computes the quantiles of the bivariate tail copula,
\eqn{r^{-1}_{2|1}(w_2|u_1=x_1;\theta)} for \eqn{0<w_2<1}.
The available tail copula families are:
\itemize{
\item 1=Husler-Reiss
\item 2=Negative logistic
\item 3=Logistic
\item 4=Dirichlet
For details on the explicit form of tail copula densities, refer to Appendix F.1 in Kiriliouk, A., Lee, J., & Segers, J. (2023).
}
}
\references{
Kiriliouk, A., Lee, J., & Segers, J. (2023). X-Vine Models for Multivariate Extremes. arXiv preprint arXiv:2312.15205.
}
