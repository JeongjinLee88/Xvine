# Internal function from the VineCopula package by Thomas Nagler and Ulf Schepsmeier and Jakob Stoeber and Eike Christian Brechmann and Benedikt Graeler and Tobias Erhardt
set_weight <- function (x, E) 
{
  is.edge <- (x[1] == E$nums[, 1]) & (x[2] == E$nums[, 2])
  if (!any(is.edge)) 
    Inf
  else -E$weights[which(is.edge)]
}
