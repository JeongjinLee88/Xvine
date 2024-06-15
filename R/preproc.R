# Internal function from the VineCopula package by Thomas Nagler and Ulf Schepsmeier and Jakob Stoeber and Eike Christian Brechmann and Benedikt Graeler and Tobias Erhardt
preproc <- function (args, ..., na.txt = NULL) 
{
  args$na.txt <- na.txt
  args$n <- length(args$u1)
  tasks <- list(...)
  for (i in seq_along(tasks)) {
    stopifnot(is.function(tasks[[i]]))
    args <- tasks[[i]](args)
  }
  args
}