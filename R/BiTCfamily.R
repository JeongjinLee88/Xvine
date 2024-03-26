#' Quantile function of a bivariate tail copula
#'
#' @description
#' Conditioning on the second argument of the bivariate tail copula being less than 1,
#' \eqn{x_1 < 1}, the function `InvTC()` computes the quantiles of the bivariate tail copula,
#'  \eqn{r^{-1}_{2|1}(w_2|u_1=x_1;\theta)} for \eqn{0<w_2<1}.
#' The available tail copula families are:
#' * 1=Husler-Reiss
#' * 2=Negative logistic
#' * 3=Logistic
#' * 4=Dirichlet
#' For details on the explicit form of tail copula densities, refer to Appendix F.1 in Kiriliouk, A., Lee, J., & Segers, J. (2023).
#'  
#' @param w2 A numeric vector of the first argument for the bivariate tail copula.
#' @param u1 A numeric vector of the second argument conditioned on \eqn{u_1 < 1}
#'  for the bivariate tail copula.
#' @param par A numeric of parameter values.
#' @param family A numeric vector indicating the list of tail copula families.
#'
#' @return A numeric vector of quantile values.
#' @export
#'
#' @references Kiriliouk, A., Lee, J., & Segers, J. (2023). X-Vine Models for Multivariate Extremes. arXiv preprint arXiv:2312.15205.
InvTC <- function(w2,u1,par,family){
  ##  HR tail copula density
  if(family==1){
    InvExp1.2=u1*exp(stats::qnorm(w2)*sqrt(par)+0.5*par)
  }
  ##  Negative logistic model
  if(family==2){
    if(par < 1e-06){
      InvExp1.2=0
    }else{
      InvExp1.2=u1*(w2^(-(par/(par+1)))-1)^(-1/par)  
    }
  }
  ##  Logistic model
  if(family==3){
    if(par==1){
      InvExp1.2=0
    }else{
      InvExp1.2=u1*((1-w2)^(par/(1-par))-1)^(1/par)  
    }
  }
  ##  Dirichlet model
  if(family==4){
    #if(par<1e-04){
    #  InvExp1.2=1e+13
    #}else{
      InvExp1.2=u1*(stats::qbeta(w2,par+1,par)/(1-stats::qbeta(w2,par+1,par)))   
    #}
  }
  return(InvExp1.2)
}


#' Bivariate tail copula with the second argument being conditioned
#' 
#' @description
#' `CondTC()` evaluates the univariate distribution of the bivariate tail copula
#'  given one of its arguments is less than 1, \eqn{\Lambda_{1|2}(u_1=x_1|x_2;\theta)}
#'  for \eqn{0< x_2 <1}.
#' 
#' @param x1 A numeric vector of the first argument for the bivariate tail copula.
#' @param x2 A numeric vector of the second argument for the bivariate tail copula.
#' @param par A numeric value of parameter values.
#' @param family A numeric vector indicating the list of bivariate tail copulas.
#' The available tail copula families are:
#' * 1=Husler-Reiss
#' * 2=Negative logistic
#' * 3=Logistic
#' * 4=Dirichlet
#' 
#' @return A numeric vector of the bivariate tail copula with a fixed argument.
#' @export
#'
#' @references Kiriliouk, A., Lee, J., & Segers, J. (2023). X-Vine Models for Multivariate Extremes. arXiv preprint arXiv:2312.15205.
CondTC <- function(x1,x2,par,family){
  ##  Log-Gaussian dist (HR)
  if(family==1){
    CondExp1.2=stats::pnorm(q = log(x1/x2),mean = par/2,sd = sqrt(par))
  }
  ##  Negative logistic
  if(family==2){
    if(par < 1e-05){
      CondExp1.2=0
    }else{
      CondExp1.2=((x1/x2)^(-par)+1)^(-1/par - 1)  
    }
  }
  ##  Logistic model
  if(family==3){
    if(par==1){
      CondExp1.2=0
    }else{
      CondExp1.2=1-((x1/x2)^par+1)^(1/par-1)  
    }
  }
  ##  Dirichlet model
  if(family==4){
    CondExp1.2=stats::pbeta(x1/(x1+x2),par+1,par)
  }
  return(CondExp1.2)
}
