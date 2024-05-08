XVineFlightGraph <- function(XVineOut, TreeDimension, IATAS, Color4Edges)
{
  VarNum=nrow(XVineOut$MST[[1]]$E$nums)+1
  if(is.null(TreeDimension)){
    TreeSeq=length(XVineOut$MST)  
  }else{
    TreeSeq=TreeDimension  
  }
  
  EdgeNum4Ind=rep(NA,TreeDimension)
  EdgeInfo=do.call(rbind,XVineOut$MST[[1]]$E$conditionedSet)
  EdgeNum4Ind[1]=nrow(EdgeInfo)
  for(i in 2:TreeSeq){
    b=do.call(rbind,XVineOut$VineTree[[i]]$E$conditionedSet)
    Pair4Ind=XVineOut$VineTree[[i]]$E$Copula.type
    if(identical(which(Pair4Ind==0),integer(0))){
      b=b
    }else{
      b=b[-which(Pair4Ind==0),]
    }
    EdgeNum4Ind[i]=nrow(b)
    EdgeInfo=rbind(EdgeInfo,b)
  }
  
  
  if(is.null(Color4Edges)){
    Color4Edges=c(rep("royalblue1",EdgeNum4Ind[1]),rep("seagreen",EdgeNum4Ind[2]),rep("lightsalmon",EdgeNum4Ind[3]),rep("orchid1",EdgeNum4Ind[4]),rep("sienna",EdgeNum4Ind[5]),rep("brown1",EdgeNum4Ind[6]),rep("navajowhite",EdgeNum4Ind[7])) 
  }
  
  ##  Plot X-vine graphs
  FlightSuperImposed=plotFlights_adj(
    IATAS,
    graph = EdgeInfo,
    edgeColors = Color4Edges,
    xyRatio = 1,
    clipMap = 1,
    edgeAlpha = 0.4
  )
  return(FlightSuperImposed)
}

