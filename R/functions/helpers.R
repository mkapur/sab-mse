## some misc helper functions from stockassessment package
setSeq <- function (min, max) 
{
  if (min == max) {
    ret <- 1
  }
  else {
    ret <- c(1:(max - min), max - min)
  }
  return(ret)
}


setS <- function (x) 
{
  setSeq(1, length(x))
}

getLowerBounds <- function (parameters) 
{
  list(sigmaObsParUS = rep(-10, length(parameters$sigmaObsParUS)))
}
