## based on stockassessment::decfon

modConfig <- function (dat) 
{
  fleetTypes <- dat$fleetTypes
  ages <- do.call(rbind, tapply(dat$aux[, 3], INDEX = dat$aux[, 
                                                              2], FUN = range))
  ages[fleetTypes %in% c(3, 5), ] <- NA
  minAge <- min(ages, na.rm = TRUE)
  maxAge <- max(ages, na.rm = TRUE)
  ages[is.na(ages)] <- minAge
  nAges <- maxAge - minAge + 1
  nFleets <- nrow(ages)
  ret <- list()
  ret$minAge <- minAge
  ret$maxAge <- maxAge
  ret$maxAgePlusGroup <- as.integer(ages[, 2] == max(ages[, 
                                                          2], na.rm = TRUE))
  x <- matrix(0, nrow = nFleets, ncol = nAges)
  lastMax <- 0
  for (i in 1:nrow(x)) {
    if (fleetTypes[i] == 0) {
      aa <- ages[i, 1]:ages[i, 2]
      aa <- aa[tapply(dat$logobs[dat$aux[, 2] == i], INDEX = dat$aux[, 
                                                                     3][dat$aux[, 2] == i], function(x) !all(is.na(x)))]
      x[i, aa - minAge + 1] <- setS(aa) + lastMax
      lastMax <- max(x)
    }
  }
  ret$keyLogFsta <- x - 1
  ret$corFlag <- 2
  x <- matrix(0, nrow = nFleets, ncol = nAges)
  lastMax <- 0
  for (i in 1:nrow(x)) {
    if (fleetTypes[i] %in% c(1, 2, 3)) {
      x[i, (ages[i, 1] - minAge + 1):(ages[i, 2] - minAge + 
                                        1)] <- setSeq(ages[i, 1], ages[i, 2]) + lastMax
      lastMax <- max(x)
    }
  }
  ret$keyLogFpar <- x - 1
  ret$keyQpow <- matrix(-1, nrow = nFleets, ncol = nAges)
  x <- matrix(0, nrow = nFleets, ncol = nAges)
  lastMax <- 0
  for (i in 1:nrow(x)) {
    if (fleetTypes[i] == 0) {
      x[i, (ages[i, 1] - minAge + 1):(ages[i, 2] - minAge + 
                                        1)] <- lastMax + 1
      lastMax <- max(x)
    }
  }
  ret$keyVarF <- x - 1
  ret$keyVarLogN <- c(1, rep(2, nAges - 1)) - 1
  x <- matrix(0, nrow = nFleets, ncol = nAges)
  lastMax <- 0
  for (i in 1:nrow(x)) {
    if (fleetTypes[i] %in% c(0, 1, 2, 3)) {
      x[i, (ages[i, 1] - minAge + 1):(ages[i, 2] - minAge + 
                                        1)] <- lastMax + 1
      lastMax <- max(x)
    }
  }
  ret$keyVarObs <- x - 1
  ret$obsCorStruct <- factor(rep("ID", nFleets), levels = c("ID", 
                                                            "AR", "US"))
  ret$keyCorObs <- matrix(-1, nrow = nFleets, ncol = nAges - 
                            1)
  colnames(ret$keyCorObs) <- paste(minAge:(maxAge - 1), (minAge + 
                                                           1):maxAge, sep = "-")
  for (i in 1:nrow(x)) {
    if (ages[i, 1] < ages[i, 2]) {
      ret$keyCorObs[i, (ages[i, 1] - minAge + 1):(ages[i, 
                                                       2] - minAge)] <- NA
    }
  }
  ret$stockRecruitmentModelCode <- 0
  ret$noScaledYears <- 0
  ret$keyScaledYears <- numeric(0)
  ret$keyParScaledYA <- array(0, c(0, 0))
  cs <- colSums(dat$catchMeanWeight)
  ii <- min(which(dat$fleetTypes == 0))
  tc <- tapply(dat$logobs[dat$aux[, 2] == ii], INDEX = dat$aux[, 
                                                               3][dat$aux[, 2] == ii], function(x) sum(x, na.rm = TRUE))
  tc <- tc * cs[names(cs) %in% names(tc)]
  pp <- tc/sum(tc)
  ret$fbarRange <- c(min(which(cumsum(pp) >= 0.25)), length(pp) - 
                       min(which(cumsum(rev(pp)) >= 0.25)) + 1) + (minAge - 
                                                                     1)
  ret$keyBiomassTreat <- ifelse(dat$fleetTypes == 3, 0, -1)
  ret$obsLikelihoodFlag <- factor(rep("LN", nFleets), levels = c("LN", 
                                                                 "ALN"))
  ret$fixVarToWeight <- 0
  ret$fracMixF <- 0
  ret$fracMixN <- 0
  ret$fracMixObs <- rep(0, nFleets)
  ret$constRecBreaks <- numeric(0)
  ret$predVarObsLink <- matrix(NA, nrow = nFleets, ncol = maxAge - 
                                 minAge + 1)
  for (i in 1:nrow(x)) {
    if (ages[i, 1] < ages[i, 2]) {
      ret$predVarObsLink[i, (ages[i, 1] - minAge + 1):(ages[i, 
                                                            2] - minAge + 1)] <- -1
    }
  }
  return(ret)
}
