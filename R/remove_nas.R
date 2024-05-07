remove_nas <- function (args) 
{
  if (!is.null(args$data)) {
    if (any(is.na(args$data))) {
      num.na <- sum(!complete.cases(args$data))
      freq.na <- round(num.na/nrow(args$data) * 100, 1)
      warning(" In ", args$call[1], ": ", num.na, " observation", 
              ifelse(num.na > 1, "s", ""), " (", freq.na, "%) contain", 
              ifelse(num.na == 1, "s", ""), " NAs.", args$na.txt, 
              call. = FALSE)
      args$na.ind <- which(!complete.cases(args$data))
      args$data <- args$data[complete.cases(args$data), 
                             , drop = FALSE]
      args$n <- nrow(args$data)
    }
  }
  else {
    if (any(is.na(args$u1 + args$u2))) {
      num.na <- sum(!complete.cases(args$u1 + args$u2))
      freq.na <- round(num.na/length(args$u1) * 100, 1)
      warning(" In ", args$call[1], ": ", num.na, " observation", 
              ifelse(num.na > 1, "s", ""), " (", freq.na, "%) contain", 
              ifelse(num.na == 1, "s", ""), " NAs.", args$na.txt, 
              call. = FALSE)
      args$na.ind <- which(is.na(args$u1 + args$u2))
      args$u1 <- args$u1[-args$na.ind]
      args$u2 <- args$u2[-args$na.ind]
    }
    else {
      args$msg <- na.ind <- NULL
    }
    args$n <- length(args$u1)
  }
  args
}