#' Converting ML estimates to dependence measures
#' 
#' @description
#' `ML2DepMtx` takes a family matrix and parameter matrix and
#' converts ML estimates for X-vine models to dependence measures via parametric relationships between
#' parameters associated with (tail) copula models and (tail) dependence measures.
#'  The first row includes the (lower) tail dependence coefficient \eqn{\chi_{L,e}=R(1,1)} for \eqn{e\in T_1}.
#'  Possible tail copula families include:
#'  Husler-Reiss, Negative logistic, logistic, and Dirichlet models. 
#'  For explicit formulas, refer to Kiriliouk, A., Lee, J., & Segers, J. (2023).
#'   and the remaining rows include the measure of Kendall's tau, \eqn{\tau_{e}}, for \eqn{e\in T_i}, \eqn{i=2,\ldots,d-1}.
#'  Any copula family with a single parameter can be considered. For more details, refer to Czado (2019, Section 3.5).
#'
#' @param FamMtx A \eqn{d\times d} upper triangular family matrix indicating the family of bivariate (tail) copula models
#' @param ParMtx A \eqn{d\times d} upper triangular parameter matrix including ML estimates for each edge
#'
#' @return A \eqn{d\times d} upper triangular matrix of dependence measures 
#' @export
#' 
#' @references Kiriliouk, A., Lee, J., & Segers, J. (2023). X-Vine Models for Multivariate Extremes. arXiv preprint arXiv:2312.15205.
#'  Czado, C. (2019). Analyzing dependent data with vine copulas. Lecture Notes in Statistics, Springer, 222.
#'
#' @examples
#' ParamM <- matrix(c(0, 1.5, 2, 2.5, 2,
#'                 0, 0, 2, 2.5, 0.7,
#'                 0, 0, 0, 0.4, -0.3,
#'                 0, 0, 0, 0, 0.1,
#'                 0, 0, 0, 0, 0),5,byrow = TRUE)
#' FamilyM <- matrix(c(0, 1, 2, 3, 4,
#'                    0, 0, 3, 4, 1,
#'                    0, 0, 0, 3, 1,
#'                    0, 0, 0, 0, 1,
#'                    0, 0, 0, 0, 0),5,byrow = TRUE)
#' ML2DepMtx(FamMtx=FamilyM,ParMtx=ParamM)               
ML2DepMtx <- function(FamMtx,ParMtx)
{
  d=dim(FamMtx)[1]
  fam1 <- FamMtx
  par <- ParMtx
  DepMtx=matrix(0,d,d)
  
  for(i in 1:(d-1)){
    if(i==1){
      for(j in (i+1):d){
        if(fam1[1,j]==1){
          DepMtx[i,j]=2-2*pnorm(sqrt(par[1,j])/2)
        }
        if(fam1[1,j]==2){
          DepMtx[i,j]=2^(-1/par[1,j])
        }
        if(fam1[1,j]==3){
          DepMtx[i,j]=2-2^(1/par[1,j])
        }
        if(fam1[1,j]==4){
          DepMtx[i,j]=integrate(Vectorize(function(x)pbeta(1/(1+x),par[1,j]+1,par[1,j])),lower = 0,upper = 1)$value
        }
      }
    }
    if(i > 1){
      for(j in (i+1):d){
        DepMtx[i,j]=VineCopula::BiCopPar2Tau(family = fam1[i,j],par = par[i,j])
      }
    }
  }
  return(DepMtx)
}
