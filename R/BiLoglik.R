#' Log-likelihood function for bivariate tail copula densities on product space
#' 
#' @description
#' `LL.BiTC` calculates a negative log-likelihood of a bivariate tail copula density \eqn{r} on product space with uniform margins.
#' As the bivariate tail copula \eqn{r(x1,x2)} with \eqn{x2<1} corresponds to the bivariate inverted-Pareto density,
#' the normalizing constant can be ignored because it is 1.
#' This function is passed to the function \code{mleBiTC()} for finding ML estimates. 
#' 
#' @param par Numeric; a parameter for bivariate tail copulas
#' @param x A \eqn{N\times 2} data matrix from a bivariate inverted-Pareto distribution
#' where the second column is less than 1
#' @param family An integer; indicates the type of bivariate tail copula families.
#' Possible tail copula families include:
#' 1=Husler-Reiss
#' 2=Negative Logistic
#' 3=Logistic
#' 4=Dirichlet
#'
#' @return Numeric; a negative log-likelihood value
#' @export
LL.BiTC <- function(par, x, family){
  n=nrow(x)
  x_c=x[,1];x_d=x[,2]
  ##  log-likelihood for HR
  if(family==1){
    ll=sum(log(x_c^(-1)*(1/sqrt(2*pi*par))*exp(-(2*par)^(-1)*(log(x_c/x_d)-par/2)^2)))
    #ll=sum(log(x_c^(-1)*(1/sqrt(2*pi*par))*exp(-(2*par)^(-1)*(log(x_c/x_d)-par/2)^2)/(2-2*pnorm(sqrt(par)/2))))
    #ll=-(n/2)*log(par)-(1/(2*par))*sum((log(x_c/x_d)-par/2)^2) # Note that there is no the extremal coefficient
  }
  ##  log-likelihood for Negative logistic
  if(family==2){
    #ll=sum(log((1+par)*(x_c*x_d)^(-par-1)*(x_c^(-par)+x_d^(-par))^(-1/par-2)))
    ll=n*log(par+1)-(par+1)*sum(log(x_c*x_d))-(1/par+2)*sum(log(x_c^(-par)+x_d^(-par)))
  }
  ##  log-likelihood for logistic
  if(family==3){
    ll=n*log(par-1)+(par-1)*sum(log(x_c*x_d))+(1/par-2)*sum(log(x_c^(par)+x_d^(par)))
  }
  ##  log-likelihood for Dirichlet
  if(family==4){
    B=gamma(par+1)*gamma(par)/gamma(2*par+1)
    ll=-n*log(B)-(2*par+1)*sum(log(x_c+x_d))+par*sum(log(x_c*x_d))
  }
  if (is.finite(ll)) {
    return(-ll) #returns the negative log-likelihood value
  }else {
    return(-10^305)
  }
}

#' Log-likelihood function for bivariate inverted-Pareto density
#'
#' @description
#' `LL.BiInvPa` calculates a negative log-likelihood of a bivariate inverted-Pareto density.
#' For details on the inverted-Pareto distribution, refer to Kiriliouk, A., Lee, J., & Segers, J. (2023).
#' 
#' 
#' @param par Numeric; a parameter for bivariate inverted-Pareto densities
#' @param x A \eqn{N\times 2} data matrix from a bivariate inverted-Pareto distribution
#' @param family An integer; indicates the type of bivariate tail copula families.
#' Possible tail copula families include:
#' 1=Husler-Reiss
#' 2=Negative Logistic
#' 3=Logistic
#' 4=Dirichlet
#'
#' @return Numeric; a negative log-likelihood value
#' @export
#' 
#' @references Kiriliouk, A., Lee, J., & Segers, J. (2023). X-Vine Models for Multivariate Extremes. arXiv preprint arXiv:2312.15205.
LL.BiInvPa <- function(par, x, family){
  n=nrow(x)
  x_c=x[,1];x_d=x[,2]
  ##  log-likelihood for HR
  if(family==1){
    ll=sum(log(x_c^(-1)*(1/sqrt(2*pi*par))*exp(-(2*par)^(-1)*(log(x_c/x_d)-par/2)^2)/(2-2*pnorm(sqrt(par)/2))))
    #ll=graphicalExtremes:::loglik_HR(data = 1/x,Gamma = par,cens = F)[1] # switch to Pareto scale / check if they are same
  }
  ##  log-likelihood for Negative logistic
  if(family==2){
    ll=n*log(par+1)-(par+1)*sum(log(x_c*x_d))-(1/par+2)*sum(log(x_c^(-par)+x_d^(-par)))-n*log(2-2^(-1/par))
  }
  ##  log-likelihood for logistic
  if(family==3){
    ll=n*log(par-1)+(par-1)*sum(log(x_c*x_d))+(1/par-2)*sum(log(x_c^(par)+x_d^(par)))-n*log(2^(1/par))
  }
  ##  log-likelihood for Dirichlet
  if(family==4){
    B=gamma(par+1)*gamma(par)/gamma(2*par+1)
    ll=-n*log(B)-(2*par+1)*sum(log(x_c+x_d))+par*sum(log(x_c*x_d))
  }
  if (is.finite(ll)) {
    return(ll) #returns the loglik
  }else {
    return(10^305)
  }
}


