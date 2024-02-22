SecondArg <- function(d,Rmat){
  CondMat=matrix(rep(NA,d^2),nrow = d)
  for(i in 2:d){
    for(k in (i-1):1){
      diff=i-k
      if(k==1 & i==2){
        CondMat[k,i]=Rmat[k, i-1] == Rmat[k,i]
      }else if(k > 1){
        if(diff == 1){
          CondMat[k,i]=Rmat[k,i-1]==Rmat[k,i]
        }else if(diff > 1){
            CondMat[k,i]=Rmat[i-1,i-1]==Rmat[k,i]
        }
        }
    }
    }
  return(CondMat)
}

SecondArg2 <- function(d,Rmat){
  CondMat2=matrix(rep(NA,d^2),nrow = d)
  for(i in 3:d){
    for(k in (i-1):2){
      diff=i-k
      if(diff > 1){
            a=(i-1):1
            CondMat2[k,i]=setequal(x=Rmat[a[-(2:(i-k))],i-1],y=Rmat[k:1,i])
        }
      }
    }
  return(CondMat2)
}

SecondArg3 <- function(d, Rmat){
  CondMat3=array(rep(NA,d*d*d),dim = c(d,d,d))
    for(i in 4:d){
      for(k in (i-2):2){
        for(j in (i-2):k){
        diff=j-k
        if(diff == 0){
          CondMat3[k,j,i]=setequal(x=Rmat[k:1,k],y=Rmat[k:1,i])
        }else if (diff==1){
          a=j:1
          b=1:diff+1
          CondMat3[k,j,i]=setequal(x=Rmat[a[-b],j],y=Rmat[k:1,i])
        }else if (diff > 1){
          a=j:1
          b=1:diff+1
          CondMat3[k,j,i]=setequal(x=Rmat[a[-b],j],y=Rmat[k:1,i])
        }
      }
    }
    }
  return(CondMat3)
}

SecondArg4 <- function(d, Rmat){
  CondMat4=array(rep(NA,d*d*d),dim = c(d,d,d))
  for(i in 4:d){
    for(k in (i-2):2){
      for(j in (i-2):k){
        diff=j-k
        if(diff == 0){
          CondMat4[k,j,i]=Rmat[k:1,k][1]==Rmat[k:1,i][1]
        }else if (diff==1){
          a=j:1
          b=1:diff+1
          CondMat4[k,j,i]=Rmat[a[-b],j][1]==Rmat[k:1,i][1]
        }else if (diff > 1){
          a=j:1
          b=1:diff+1
          CondMat4[k,j,i]=Rmat[a[-b],j][1]==Rmat[k:1,i][1]
        }
      }
    }
  }
  return(CondMat4)
}

