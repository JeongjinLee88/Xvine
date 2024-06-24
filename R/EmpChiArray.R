EmpChiArray <- function (Dat, quan = NULL) 
{
  if (!is.matrix(Dat)) {
    stop("The data should be a matrix")
  }
  if (ncol(Dat) <= 1) {
    stop("The data should be a matrix with at least two columns.")
  }
  if (!is.null(quan)) {
    # switch to uniform scale
    Dat <- ParetoTransRank(data = Dat, u_quan = quan, scaleType = "U")
  }
  d <- ncol(Dat)
  ChiArray=array(NA,dim = c(d,d,d))
  for(k in 1:d){
    for(i in 1:d){
      for(j in 1:d){
        Z=Dat[,c(i,j,k)] # 123
        Z12.3=Z[Z[,3]<1,] #(X1,X2,X3) given X3<1
        Z13.2=Z[Z[,2]<1,c(1,3,2)] #(X1,X3,X2) given X2<1
        Z23.1=Z[Z[,1]<1,c(2,3,1)] #(X2,X3,X1) given X1<1
        
        chi1 <- sum(apply(Z12.3[,1:2],1,function(x)all(x<1)))/nrow(Z12.3) # X1<1,X2<1 | X3<1
        chi2 <- sum(apply(Z13.2[,1:2],1,function(x)all(x<1)))/nrow(Z13.2) # X1<1,X3<1 | X2<1
        chi3 <- sum(apply(Z23.1[,1:2],1,function(x)all(x<1)))/nrow(Z23.1) # X2<1,X3<1 | X1<1
        chi <- mean(c(chi1,chi2,chi3))
        ChiArray[i,j,k]=chi 
      }
    }
  }
  return(ChiArray)
}
