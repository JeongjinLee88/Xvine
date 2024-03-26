#' Redefining the conditional distribution function of pair-copulas
#' 
#' @description
#' `PseudoCop` recalculates the conditional distribution function of pair-copulas by conditioning on all variables
#' in the conditioning set being less than 1.
#' This calculation is based on recursive approaches tree by tree sequentially from \eqn{T_3} to \eqn{T_{d-1}}.
#' For more details, refer to Kiriliouk, A., Lee, J., & Segers, J. (2023)
#' 
#' @param j Integer; the number of conditioning variables in the conditioning set \eqn{D(e)}
#' @param cind Integer; stands for the column index corresponding to the \eqn{i}th column. 
#' Given the \eqn{i}th column, subset the data given all variables in \eqn{D(e)} are less than 1.
#' @param data An \eqn{N\times d} data matrix of either pseudo-observations or random samples from the inverted-Pareto distribution
#' @param StrMtx A \eqn{d \times d} specified structure matrix
#' @param fam A \eqn{d \times d} specified family matrix
#' @param MaxMtr A \eqn{d\times d} max-matrix
#' @param Par A \eqn{d\times d} recursively estimated parameter matrix from previous trees
#'
#' @references Kiriliouk, A., Lee, J., & Segers, J. (2023). X-Vine Models for Multivariate Extremes. arXiv preprint arXiv:2312.15205.
#' @return A list of two numeric vectors, including two conditional distribution functions of pair-copulas:
#'  * VdirC
#'  * VindirC
#' @export
PseudoCop <- function(j, cind, data, StrMtx, fam, MaxMtr, Par){
  d <- ncol(data)
  # if j=2 (the number of conditioning variables in D(e)), then
  # k-> l=1,2 tree levels
  # i-> s=5,4,3,2 / 5,4,3 / 5,4 each column
  VdirC <- list()
  VindirC <- list()
  for(l in 1:j){
    if(l==1){
      for(s in d:(l+1)){ 
        ind=StrMtx[1:j,cind] #does not depend on "l"
        #CondData <- data[data[,ind[1]]<1 & data[,ind[2]]<1 ,]
        CondData <- data[apply(data[,c(ind)],1,function(x)all(x<1)),]
        m <- MaxMtr[l, s]
        Z=CondData[,c(s,m)] # 54
        Z1.2=Z[Z[,2]<1,] #lambda5.4 given X4<1
        Z2.1=Z[Z[,1]<1,c(2,1)] #lambda4.5 given X5<1
        
        direct <- CondTC(x1 = Z1.2[,1],x2 = Z1.2[,2],par = Par[l,s],family = fam[l,s])
        indirect <- CondTC(x1 = Z2.1[,1],x2 = Z2.1[,2],par = Par[l,s],family = fam[l,s])
        VdirC[[s]]=direct
        VindirC[[s]]=indirect 
      }
    }
    if(l > 1){
      for(s in d:(l+1)){ 
        m <- MaxMtr[l, s]
        zr1 <- VdirC[[s]] 
        zr2 <- if (m == StrMtx[l, s]) {
          VdirC[[m]]
        }else {
          VindirC[[m]] 
        }
        Z=cbind(zr1,zr2)
        #direct <- CondCop(u1 = zr1,u2 = zr2,par = Par[l,s],family = fam[l,s])
        #indirect <- CondCop(u1 = zr2,u2 = zr1,par = Par[l,s],family = fam[l,s])
        direct <- BiCopHfunc2(u1 = zr1,u2 = zr2,par = Par[l,s],family = fam[l,s])
        indirect <- BiCopHfunc2(u1 = zr2,u2 = zr1,par = Par[l,s],family = fam[l,s])
        VdirC[[s]]=direct
        VindirC[[s]]=indirect
      }  
    }
  }
  return(list("Vdir"=VdirC,"Vindir"=VindirC))
}

