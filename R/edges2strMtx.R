edgesToM <- function(edges, format = F){ 
  # add check for edges format?
  if(format == T){ 
    edges <- append(list(t(edges[[1]][[1]])), lapply(edges[-1], function(i) rbind(t(i[[1]]),t(i[[2]]))))
    names(edges) <- NULL
  }
  
  l <- length(edges) 
  d <- length(unique(c(edges[[1]])))
  M <- matrix(0, nrow = d, ncol = d)
  for(i in d:2){
    # we have to pick an index so that after removing all nodes with this value, the result is still a vine
    inds <- c(edges[[l]][1:2,])
    temp <- unique(as.numeric(names(which(table(inds) == 1))))
    if(i > 2){
      M[i,i] <- max(temp[which(!(temp %in% unique(c(edges[[l]][3:(l+1),]))))])
    } else{
      M[i,i] <- max(temp)
    }
    indices <- sapply(edges, function(mat) which(mat[1:2,,drop=F] == M[i,i], arr.ind = T)[1,2])
    temp <- unlist(sapply(c(1:l), function(j) edges[[j]][1:2,indices[j]]))
    diagn <- temp[temp!=M[i,i]]
    if(length(diagn) != (i-1)){
      diagn <- c(diagn,rep(0,(i-1)-length(diagn)))
    }
    M[1:(i-1),i] <- diagn
    edges <- lapply(edges, function(mat) {  # removes all nodes which contain the chosen index: 
      mat <- mat[,apply(mat, 2, function(col) all(col[1:2] != M[i,i])), drop = F]
    })
    if(dim(edges[[length(edges)]])[2] == 0){
      edges <- edges[-length(edges)]
      l <- l-1
    }
  }
  M[1,1] <- setdiff(c(1:d), diag(M)[2:d])
  return(M)
}

#################################
########### Standard example
##################################
M <- cbind(c(1,0,0,0,0),c(1,2,0,0,0),c(2,1,3,0,0),c(2,3,1,4,0),c(4,2,3,1,5))
(edges <- XVineEdges(M)$edges)
edgesToM(edges) #ok! 

(edges <- XVineEdges(M, format = T)$edges)
edgesToM(edges, format = T) #ok! 

# Same vine structure but different initial M (since it is not unique) 
M2 <- cbind(c(2,0,0,0,0),c(2,4,0,0,0),c(4,2,5,0,0),c(2,4,5,3,0),c(2,3,4,5,1))
(edges <- XVineEdges(M2)$edges)
edgesToM(edges) #ok!

## Same example but truncated after T2
M <- cbind(c(1,0,0,0,0),c(1,2,0,0,0),c(2,1,3,0,0),c(2,3,0,4,0),c(4,2,0,0,5))
(edges <- XVineEdges(M)$edges)
edgesToM(edges) #ok!
(edges <- XVineEdges(M, format = T)$edges)
edgesToM(edges, format = T)

# Same vine structure but different initial M (since it is not unique) 
M2 <- cbind(c(2,0,0,0,0),c(2,4,0,0,0),c(4,2,5,0,0),c(2,4,0,3,0),c(2,3,0,0,1))
(edges <- XVineEdges(M2)$edges)
edgesToM(edges) #ok!

# With Mmod and Mpar (models & parameter values randomly chosen)
# Do we need to 'format' these as well?
Mmod <- cbind(rep(0,5),c("neglog", rep(0,4)), c("hr", "clayton", 0, 0, 0), 
              c("log", "gumbel", 0, 0, 0), c("hr", "frank", 0, 0, 0))
Mpar <- cbind(rep(0,5),c(0.5, rep(0,4)), c(0.7, 0.3, 0, 0, 0), 
              c(0.9,0.1,0, 0, 0), c(0.55,0.25,0,0, 0))
XVineEdges(M, Mmod, Mpar)


#################################
########### Example 5.11 on page 112 of the vine copula book
##################################

M <- cbind(c(4,0,0,0,0,0), c(4,3,0,0,0,0), c(3,4,1,0,0,0), c(1,3,4,2,0,0), 
           c(1,3,4,2,5,0), c(5,1,3,4,2,6))
(edges <- XVineEdges(M)$edges)
M2 <- edgesToM(edges)  #gives a different M matrix than the initial one, but is correct! 
(edges2 <- XVineEdges(M2)$edges) 

# Same vine structure but different initial M (since it is not unique) 
M3 <- cbind(c(5,0,0,0,0,0), c(5,1,0,0,0,0), c(5,1,6,0,0,0), c(1,5,6,3,0,0), 
            c(3,1,5,6,4,0), c(1,3,4,5,6,2))
(edges <- XVineEdges(M)$edges)
M4 <- edgesToM(edges)  #gives a different M matrix than the initial one, but is correct! 
(edges4 <- XVineEdges(M4)$edges) 


## Same example but truncated after T1
M <- cbind(c(4,0,0,0,0,0), c(4,3,0,0,0,0), c(3,0,1,0,0,0), c(1,0,0,2,0,0), 
           c(1,0,0,0,5,0), c(5,0,0,0,0,6))
(edges <- XVineEdges(M)$edges)
(edges <- XVineEdges(M, format = T)$edges)

## Same example but truncated after T3
M <- cbind(c(4,0,0,0,0,0), c(4,3,0,0,0,0), c(3,4,1,0,0,0), c(1,3,4,2,0,0), 
           c(1,3,4,0,5,0), c(5,1,3,0,0,6))
(edges <- XVineEdges(M)$edges)
M2 <- edgesToM(edges)  #gives a different M matrix than the initial one, but is correct! 
(edges2 <- XVineEdges(M2)$edges) # of course, the order of the columns in each element of 'edges' is different


## Same example but truncated after T5
M <- cbind(c(4,0,0,0,0,0), c(4,3,0,0,0,0), c(3,4,1,0,0,0), c(1,3,4,2,0,0), 
           c(1,3,4,2,5,0), c(5,1,3,4,0,6))
(edges <- XVineEdges(M)$edges)
M2 <- edgesToM(edges)  #gives a different M matrix than the initial one, but is correct! 
(edges2 <- XVineEdges(M2)$edges) # of course, the order of the columns in each element of 'edges' is different


#################################
########### Jeongjin's example in d = 6
##################################
M <- cbind(c(1,0,0,0,0,0), c(1,4,0,0,0,0), c(4,1,6,0,0,0), c(4,1,6,3,0,0), 
           c(6,4,1,3,5,0), c(3,4,1,6,5,2))
(edges <- XVineEdges(M, format = T)$edges) #ok 
(M2 <- edgesToM(edges, format = T))  #gives a different M matrix than the initial one, but is correct! 
(edges2 <- XVineEdges(M2, format = T)$edges) 

#truncation after level 3
M <- cbind(c(1,0,0,0,0,0), c(1,4,0,0,0,0), c(4,1,6,0,0,0), c(4,1,6,3,0,0), 
           c(6,4,1,0,5,0), c(3,4,1,0,0,2))
(edges <- XVineEdges(M)$edges) #ok 
(M2 <- edgesToM(edges)) #ok
(edges2 <- XVineEdges(M2)$edges) 




