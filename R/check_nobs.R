# Internal function from the VineCopula package by Thomas Nagler and Ulf Schepsmeier and Jakob Stoeber and Eike Christian Brechmann and Benedikt Graeler and Tobias Erhardt
check_nobs <- function (args) 
{
  if (!is.null(args$data)) {
    if (nrow(args$data) < 2) 
      stop("\n In ", args$call[1], ": ", "Number of observations has to be at least 2.", 
           call. = FALSE)
  }
  else {
    if (length(args$u1) < 2) 
      stop("\n In ", args$call[1], ": ", "Number of observations has to be at least 2.", 
           call. = FALSE)
  }
  args
}