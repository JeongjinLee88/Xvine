#' Range of parameters for tail copula families
#' 
#' @description
#' `ParRangeTC` specifies a proper range of parameters for tail copula families.
#'  This range specification helps the convergence of the `Optimize` command.
#' 
#' @param family Integer; indicates the type of bivariate tail copula densities.
#' Possible tail copula families include:
#' 1=Husler-Reiss
#' 2=Negative Logistic
#' 3=Logistic
#' 4=Dirichlet
#'
#' @return A numeric vector; specifies the range of parameter for the corresponding tail copula density
#' @export
#'
ParRangeTC <- function(family){
  ##  Husler-Reiss
  if(family==1){
    range=c(0.01, 20)
  }
  ##  Negative logistic
  if(family==2){
    range=c(0.01, 20)
  }
  ##  Logistic
  if(family==3){
    range=c(1.01, 20)
  }
  ##  Dirichlet
  if(family==4){
    range=c(0.01, 20)
  }
  return(range)
}


