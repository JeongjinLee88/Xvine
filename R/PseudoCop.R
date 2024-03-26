PseudoCop <- function(j, cind, data, StrMtx, fam, MaxMtr, Par){
  d <- ncol(data)
  #M <- XVS$xmat[,,1]
  #M <- ReorderStrMtx(XVS$xmat[,,1])
  #fam1 <- XVS$fmat[,,1]
  #MaxMat <- createMaxMtx(Matrix = M)
  # For example, if j=2 the number of conditioning variables in D(e), then
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

# Arguments
# j: # of conditioning variables in D(e)
# cind: Integer: 'cind'(column index) indicates the 'i'th column. Given the ith column, subset the data given all variables in D(e) are less than 1
# XVS: An XVineSpec object that contains all information on the X-vine specification up to permutations
# data: A Nxd data matrix
# Par: A dxd estimated parameter matrix

# Value
# Returns redefined pseudo data either in V1 or in V2

#PseudoCop <- function(j,cind,d,Stmtx,MaxMat,data,Params,Fam){
  # For example, if j=2 the number of conditioning variables in D(e), then
  # k-> l=1,2 tree levels (the 'l'th row)
  # i-> s=5,4,3,2 / 5,4,3 / 5,4 each column
#  VdirC <- list()
#  VindirC <- list()
#  for(l in 1:j){
#    if(l==1){
#      for(s in d:(l+1)){ 
#        ind=Stmtx[1:j,cind] #does not depend on "l"
        #CondData <- data[data[,ind[1]]<1 & data[,ind[2]]<1 ,]
#        CondData <- data[apply(data[,c(ind)],1,function(x)all(x<1)),]
#        m <- MaxMat[l, s]
#        Z=CondData[,c(s,m)] # 54
#        Z1.2=Z[Z[,2]<1,] #lambda5.4 given X4<1
#        Z2.1=Z[Z[,1]<1,c(2,1)] #lambda4.5 given X5<1
        
#        direct <- CondExp(x1 = Z1.2[,1],x2 = Z1.2[,2],par = Params[l,s],family = Fam[l,s])
#        indirect <- CondExp(x1 = Z2.1[,1],x2 = Z2.1[,2],par = Params[l,s],family = Fam[l,s])
#        VdirC[[s]]=direct
#        VindirC[[s]]=indirect 
#      }
#    }
#    if(l > 1){
#      for(s in d:(l+1)){ 
#        m <- MaxMat[l, s]
#        zr1 <- VdirC[[s]] 
#        zr2 <- if (m == Stmtx[l, s]) {
#          VdirC[[m]]
#       }else {
#          VindirC[[m]] 
#        }
#        Z=cbind(zr1,zr2)
#        direct <- CondCop(u1 = zr1,u2 = zr2,par = Params[l,s],family = Fam[l,s])
#        indirect <- CondCop(u1 = zr2,u2 = zr1,par = Params[l,s],family = Fam[l,s])
#        VdirC[[s]]=direct
#        VindirC[[s]]=indirect
#      }  
#    }
#  }
#  return(list("Vdir"=VdirC,"Vindir"=VindirC))
#}

