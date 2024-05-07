check_data <- function (args) 
{
  if (is.symbol(args$data)) 
    stop("\n In ", args$call[1], ": ", "Argument data is missing.", 
         call. = FALSE)
  if (is.null(args$data)) 
    stop("\n In ", args$call[1], ": ", "Argument data is missing.", 
         call. = FALSE)
  if (is.vector(args$data)) {
    args$data <- t(as.matrix(args$data))
  }
  else {
    args$data <- as.matrix(args$data)
  }
  if (!is.numeric(args$data) & !all(is.na(args$data))) 
    stop("\n In ", args$call[1], ": ", "Data have to be numeric.", 
         call. = FALSE)
  if (ncol(args$data) < 2) 
    stop("\n In ", args$call[1], ": ", "Dimension has to be at least 2.", 
         call. = FALSE)
  if (is.null(colnames(args$data))) 
    colnames(args$data) <- paste("V", seq.int(ncol(args$data)), 
                                 sep = "")
  args$n <- nrow(args$data)
  args$d <- ncol(args$data)
  args
}