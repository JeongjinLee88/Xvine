#' Multivariate Pareto samples from a 50-dimensional X-vine model with a C-vine structure
#'
#' Given the X-vine specification, the sample data is drawn from the function `Paretosim()`.
#' The first tree consists of the Husler-Reiss models and negative logistic models with randomly assigned parameter
#' values \eqn{\theta_{e}\in[1,2]} for \eqn{e\in E_1}.
#' Subsequent trees contain bivariate Gaussian copulas with partial correlations \eqn{\rho_{e}=1.1-0.1j} for \eqn{e\in E_j} and \eqn{j\in\{2,\ldots,9\}},
#' and \eqn{\rho_{e}=0.1} for \eqn{e\in E_j} with \eqn{j\ge 10}. 
#'
#' @format A matrix with 5,000 rows and 50 columns with Pareto scale.
#' @docType data
#' @name XVinePa50dim
"XVinePa50dim"


