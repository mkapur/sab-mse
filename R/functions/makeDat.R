## based on stockassessment::setup.sam.data

makeDat <-
  function (fleets = NULL,
            surveys = NULL,
            # residual.fleet = NULL, ## we dont sep recreational
            prop.mature = NULL, ## later index by k
            stock.mean.weight = NULL,
            catch.mean.weight = NULL,
            dis.mean.weight = NULL,
            land.mean.weight = NULL,
            natural.mortality = NULL,
            prop.f = NULL,
            prop.m = NULL,
            land.frac = NULL,
            recapture = NULL)

{
  fleet.idx <- 0
  type <- NULL
  time <- NULL
  name <- NULL
  corList <- list()
  ## columns are years, 1 row per fleet
  idxCor <- matrix(NA, nrow = length(fleets) + length(surveys) + 
                     1, ncol = nrow(natural.mortality))
  colnames(idxCor) <- rownames(natural.mortality)
  dat <- data.frame(year = NA, fleet = NA, age = NA, aux = NA)
  weight <- NULL
  doone <- function(m) { ## single fleet per input dataset
    if('Yr' %in% colnames(m)  | 'Year'  | colnames(m)) year <- m$
    year <- rownames(m)[row(m)] ## rowNAMES are years
    fleet.idx <<- fleet.idx + 1 ## updates for nesting as we move along
    fleet <- rep(fleet.idx, length(year))
    age <- as.integer(colnames(m)[col(m)]) ## assumes columns are age bins
    aux <- as.vector(m) ## vectorize all contents
    dat <<- rbind(dat, data.frame(year, fleet, age, aux)) ## reshape (again)...
    if ("weight" %in% names(attributes(m))) {
      weight <<- c(weight, as.vector(attr(m, "weight")))
    }
    else {
      if ("cov" %in% names(attributes(m))) {
        weigthTmp = do.call(rbind, lapply(attr(m, "cov"), 
                                          diag))
        weight <<- c(weight, as.vector(weigthTmp))
      }
      else {
        if ("cov-weight" %in% names(attributes(m))) {
          weigthTmp = do.call(rbind, lapply(attr(m, 
                                                 "cov-weight"), diag))
          weight <<- c(weight, 1/as.vector(weigthTmp))
        }
        else {
          weight <<- c(weight, rep(NA, length(year)))
        }
      }
    }
    if ("cov" %in% names(attributes(m))) {
      attr(m, "cor") <- lapply(attr(m, "cov"), cov2cor)
    }
    if ("cov-weight" %in% names(attributes(m))) {
      attr(m, "cor") <- lapply(attr(m, "cov-weight"), 
                               cov2cor)
    }
    if ("cor" %in% names(attributes(m))) {
      thisCorList <- attr(m, "cor")
      whichCorOK <- which(unlist(lapply(thisCorList, function(x) !any(is.na(x)))))
      thisCorList <- thisCorList[whichCorOK]
      corList <<- c(corList, thisCorList)
      nextIdx <- if (all(is.na(idxCor))) {
        0
      }
      else {
        max(idxCor, na.rm = TRUE)
      }
      idxCor[fleet.idx, colnames(idxCor) %in% rownames(m)][whichCorOK] <<- nextIdx:(nextIdx + 
                                                                                      length(thisCorList) - 1)
    }
  }
  if (!is.null(residual.fleet)) {
    doone(residual.fleet) ## updates dat with this fleet info
    type <- c(type, 0)
    time <- c(time, 0)
    name <- c(name, "Residual catch")
  }
  if (!is.null(fleets)) {
    if (is.data.frame(fleets) | is.matrix(fleets)) {
      doone(fleets)
      type <- c(type, 1) ## 1 means fishery fleet
      time <- c(time, 0)
      name <- c(name, "Comm fleet")
    }
    else {
      dummy <- lapply(fleets, doone)
      type <- c(type, rep(1, length(fleets)))
      time <- c(time, rep(0, length(fleets)))
      name <- c(name, strtrim(gsub("\\s", "", names(dummy)), 
                              50))
    }
  }
  if (!is.null(surveys)) {
    if (is.data.frame(surveys) | is.matrix(surveys)) {
      doone(surveys)
      thistype <- ifelse(min(as.integer(colnames(surveys))) < 
                           (-0.5), 3, 2) ## two means normal survey (age?), three means length?
      type <- c(type, thistype)
      time <- c(time, mean(attr(surveys, "time")))
      name <- c(name, "Survey fleet")
    }
    else { ## if more than one survey loop
      dummy <- lapply(surveys, doone)
      type <- c(type, unlist(lapply(surveys, function(x) ifelse(min(as.integer(colnames(x))) < 
                                                                  (-0.5), 3, 2))))
      time <- c(time, unlist(lapply(surveys, function(x) mean(attr(x, 
                                                                   "time")))))
      name <- c(name, strtrim(gsub("\\s", "", names(dummy)), 
                              50))
    }
  }
  if (is.null(land.frac)) {
    land.frac <- matrix(1, nrow = nrow(residual.fleet), 
                        ncol = ncol(residual.fleet))
  }
  if (is.null(dis.mean.weight)) {
    dis.mean.weight <- catch.mean.weight
  }
  if (is.null(land.mean.weight)) {
    land.mean.weight <- catch.mean.weight
  }
  if (is.null(prop.f)) {
    prop.f <- matrix(0, nrow = nrow(residual.fleet), ncol = ncol(residual.fleet))
  }
  if (is.null(prop.m)) {
    prop.m <- matrix(0, nrow = nrow(residual.fleet), ncol = ncol(residual.fleet))
  }
  dat$aux[which(dat$aux <= 0)] <- NA
  dat <- dat[!is.na(dat$year), ]
  if (!is.null(recapture)) {
    tag <- data.frame(year = recapture$ReleaseY)
    fleet.idx <- fleet.idx + 1
    tag$fleet <- fleet.idx
    tag$age <- recapture$ReleaseY - recapture$Yearclass
    tag$aux <- exp(recapture$r)
    tag <- cbind(tag, recapture[, c("RecaptureY", "Yearclass", 
                                    "Nscan", "R", "Type")])
    dat[names(tag)[!names(tag) %in% names(dat)]] <- NA
    dat <- rbind(dat, tag)
    weight <- c(weight, rep(NA, nrow(tag)))
    type <- c(type, 5)
    time <- c(time, 0)
    name <- c(name, "Recaptures")
  }
  dat <- dat[complete.cases(dat[, 1:3]), ]
  o <- order(as.numeric(dat$year), as.numeric(dat$fleet), 
             as.numeric(dat$age))
  attr(dat, "type") <- type
  names(time) <- NULL
  attr(dat, "time") <- time
  names(name) <- NULL
  attr(dat, "name") <- name
  dat <- dat[o, ]
  weight <- weight[o]
  newyear <- min(as.numeric(dat$year)):max(as.numeric(dat$year))
  newfleet <- min(as.numeric(dat$fleet)):max(as.numeric(dat$fleet))
  mmfun <- function(f, y, ff) {
    idx <- which(dat$year == y & dat$fleet == f)
    ifelse(length(idx) == 0, NA, ff(idx) - 1)
  }
  idx1 <- outer(newfleet, newyear, Vectorize(mmfun, c("f", 
                                                      "y")), ff = min)
  idx2 <- outer(newfleet, newyear, Vectorize(mmfun, c("f", 
                                                      "y")), ff = max)
  attr(dat, "idx1") <- idx1
  attr(dat, "idx2") <- idx2
  attr(dat, "minAgePerFleet") <- tapply(as.integer(dat[, "age"]), 
                                        INDEX = dat[, "fleet"], FUN = min)
  attr(dat, "maxAgePerFleet") <- tapply(as.integer(dat[, "age"]), 
                                        INDEX = dat[, "fleet"], FUN = max)
  attr(dat, "year") <- newyear
  attr(dat, "nyear") <- max(as.numeric(dat$year)) - min(as.numeric(dat$year)) + 
    1
  cutY <- function(x) x[rownames(x) %in% newyear, ]
  attr(dat, "prop.mature") <- cutY(prop.mature)
  attr(dat, "stock.mean.weight") <- cutY(stock.mean.weight)
  attr(dat, "catch.mean.weight") <- cutY(catch.mean.weight)
  attr(dat, "dis.mean.weight") <- cutY(dis.mean.weight)
  attr(dat, "land.mean.weight") <- cutY(land.mean.weight)
  attr(dat, "natural.mortality") <- cutY(natural.mortality)
  attr(dat, "prop.f") <- cutY(prop.f)
  attr(dat, "prop.m") <- cutY(prop.m)
  attr(dat, "land.frac") <- cutY(land.frac)
  ret <-
    list(
      noFleets = length(attr(dat, "type")),
      fleetTypes = as.integer(attr(dat,
                                   "type")),
      sampleTimes = attr(dat, "time"),
      noYears = attr(dat,
                     "nyear"),
      years = attr(dat, "year"),
      minAgePerFleet = attr(dat,
                            "minAgePerFleet"),
      maxAgePerFleet = attr(dat, "maxAgePerFleet"),
      nobs = nrow(dat),
      idx1 = attr(dat, "idx1"),
      idx2 = attr(dat,
                  "idx2"),
      idxCor = idxCor,
      aux = do.call(cbind, lapply(dat,
                                  as.numeric))[,-4],
      logobs = log(dat[, 4]),
      weight = as.numeric(weight),
      propMat = attr(dat, "prop.mature"),
      stockMeanWeight = attr(dat,
                             "stock.mean.weight"),
      catchMeanWeight = attr(dat,
                             "catch.mean.weight"),
      natMor = attr(dat, "natural.mortality"),
      landFrac = attr(dat, "land.frac"),
      disMeanWeight = attr(dat,
                           "dis.mean.weight"),
      landMeanWeight = attr(dat, "land.mean.weight"),
      propF = attr(dat, "prop.f"),
      propM = attr(dat, "prop.m"),
      corList = corList
    )
  attr(ret, "fleetNames") <- attr(dat, "name")
  return(ret)
}
