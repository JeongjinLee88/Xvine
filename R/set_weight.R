set_weight <- function (x, E) 
{
  is.edge <- (x[1] == E$nums[, 1]) & (x[2] == E$nums[, 2])
  if (!any(is.edge)) 
    Inf
  else -E$weights[which(is.edge)]
}
