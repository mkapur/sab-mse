## based on stockassessment::clean.void.catches


cleanCatch <- function (dat, conf) 
{
  rmidx <- ((dat$aux[, 3] %in% (conf$minAge:conf$maxAge)[which(conf$keyLogFsta[1, 
                                                                               ] == (-1))]) & dat$aux[, 2] == 1)
  dat$aux <- dat$aux[!rmidx, ]
  dat$logobs <- dat$logobs[!rmidx]
  dat$weight <- dat$weight[!rmidx]
  dat$nobs <- sum(!rmidx)
  dat$minAgePerFleet <- as.integer(tapply(dat$aux[, "age"], 
                                          INDEX = dat$aux[, "fleet"], FUN = min))
  dat$maxAgePerFleet <- as.integer(tapply(dat$aux[, "age"], 
                                          INDEX = dat$aux[, "fleet"], FUN = max))
  newyear <- min(as.numeric(dat$aux[, "year"])):max(as.numeric(dat$aux[, 
                                                                       "year"]))
  newfleet <- min(as.numeric(dat$aux[, "fleet"])):max(as.numeric(dat$aux[, 
                                                                         "fleet"]))
  mmfun <- function(f, y, ff) {
    idx <- which(dat$aux[, "year"] == y & dat$aux[, "fleet"] == 
                   f)
    ifelse(length(idx) == 0, NA, ff(idx) - 1)
  }
  dat$idx1 <- outer(newfleet, newyear, Vectorize(mmfun, c("f", 
                                                          "y")), ff = min)
  dat$idx2 <- outer(newfleet, newyear, Vectorize(mmfun, c("f", 
                                                          "y")), ff = max)
  dat
}
