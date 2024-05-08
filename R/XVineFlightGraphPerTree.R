XVineFlightGraphPerTree <- function(XVineOut, IATAS, TreeLevel, Color4Edge)
{
  #VarNum=length(XVineOut$MST)+1
  if(TreeLevel==1){
    TreeColor="royalblue1"  
  }
  if(TreeLevel==2){
    TreeColor="seagreen"
  }
  if(TreeLevel==3){
    TreeColor="lightsalmon"
  }
  if(TreeLevel==4){
    TreeColor="orchid1"
  }
  if(TreeLevel==5){
    TreeColor="sienna"
  }
  if(TreeLevel==6){
    TreeColor="brown1"
  }
  if(TreeLevel==7){
    TreeColor="navajowhite"
  }
  
  EdgeInfo=do.call(rbind,XVineOut$VineTree[[TreeLevel]]$E$conditionedSet)
  Pair4Ind=XVineOut$VineTree[[TreeLevel]]$E$Copula.type
  if(identical(which(Pair4Ind==0),integer(0))){
    EdgeInfo=EdgeInfo
  }else{
    EdgeInfo=EdgeInfo[-which(Pair4Ind==0),]
  }
  EdgeN4Ind=nrow(EdgeInfo)
  
  if(is.null(Color4Edge)){
    #Color4Edge=rep(TreeColor,VarNum-TreeLevel)
    Color4Edge=rep(TreeColor,EdgeN4Ind)
  }  
  
  FlightPerTree=plotFlights_adj(
    IATAS,
    graph = EdgeInfo,
    edgeColors = Color4Edge,
    xyRatio = 1,
    clipMap = 1,
    edgeAlpha = 0.6
  )
  return(FlightPerTree)
}

