
checkDataPath <- function(path, type = c('static', 'yearly')) {

  if (!file.exists(path)) {
    return(FALSE)
  }
  ext <- tools::file_ext(path)
  if (ext != 'csv') {
    return(FALSE)
  }

  type <- match.arg(type)
  if (type == 'static') {
    staticHeader <- as.character(read.csv(file = path, header = FALSE, nrows = 1, stringsAsFactors = FALSE)[1,])
    if (!('gid' %in% staticHeader)) {
      return(FALSE)
    }
  } else {
    yearlyHeader <- as.character(read.csv(file = path, header = FALSE, nrows = 1, stringsAsFactors = FALSE)[1,])
    if (!all(c('gid', 'year') %in% yearlyHeader)) {
      return(FALSE)
    }
  }

  return(TRUE)
}

#' Set Paths to PRIO GRID csv Files
#'
#' \code{setDataPath} is used to (re)set the paths to the csv files containing the
#' PRIO GRID raw data. The path information is saved in the installation directory
#' of the package, so it is persistent across R sessions.
#'
#' It is possible to set a single path if only the static or the yearly variables are required.
#' Also, (re)setting a single path will not overwrite the other.
#'
#' @param staticPath A character string represeting a path pointing to a static PRIO GRID csv file.
#' If the default value of \code{''} is passed, no changes are made.
#' @param yearlyPath A character string represeting a path pointing to a yearly PRIO GRID csv file.
#' If the default value of \code{''} is passed, no changes are made.
#'
#' @return Void.
#'
#' @examples
#'
#' \dontrun{
#' setDataPath(staticPath = '/home/usr/data/PRIO-GRID Static Variables - 2017-12-08.csv')
#' }
#'
#' @export
setDataPath <- function(staticPath = '', yearlyPath = '') {

  ## Check contents
  if (staticPath != '') {
    staticCheck <- checkDataPath(staticPath, type = 'static')
    if (!staticCheck) {
      stop(paste(staticPath, 'is not a valid path to a static PRIO GRID csv file.'))
    }
  }
  if (yearlyPath != '') {
    yearlyCheck <- checkDataPath(yearlyPath, type = 'yearly')
    if (!yearlyCheck) {
      stop(paste(yearlyCheck, 'is not a valid path to a yearly PRIO GRID csv file.'))
    }
  }

  ## Save path object
  loc <- file.path(find.package("pgrid"), "data")
  unlink(file.path(loc, "dataPath.rda"), recursive = TRUE, force = FALSE)
  dataPath <- c(staticPath, yearlyPath)
  saveRDS(dataPath, file=file.path(loc, "dataPath.rds"))
}

isDataPathAvail <- function(type = c('static', 'yearly')) {

  type <- match.arg(type)

  loc <- file.path(find.package("pgrid"), "data")
  targetPath <- file.path(loc, "dataPath.rds")
  if (file.exists(targetPath)) {
    dataPath <- readRDS(file=targetPath)
  } else {
    stop('Paths to PRIO GRID csv files not set.')
  }

  if (type == 'static') {
    retVal <- checkDataPath(dataPath[1], type = 'static')
  } else {
    retVal <- checkDataPath(dataPath[2], type = 'yearly')
  }

  return(retVal)
}

getDataPath <- function(type = c('static', 'yearly')) {

  type <- match.arg(type)

  loc <- file.path(find.package("pgrid"), "data")
  targetPath <- file.path(loc, "dataPath.rds")
  if (file.exists(targetPath)) {
    dataPath <- readRDS(file=targetPath)
  } else {
    stop('Paths to PRIO GRID csv files not set.')
  }

  if (type == 'static') {
    staticPath <- dataPath[1]
    if (staticPath == '') {
      stop('Path to static PRIO GRID csv not set.')
    } else {
      out <- staticPath
    }
  } else {
    yearlyPath <- dataPath[2]
    if (yearlyPath == '') {
      stop('Path to yearly PRIO GRID csv not set.')
    } else {
      out <- yearlyPath
    }
  }

  return(out)
}

#' Get Available PRIO GRID Variables
#'
#' \code{getVarNames} returns a data frame of all PRIO GRID variables available to the package.
#' Which variables are available depends on the contents of the PRIO GRID csv files linked by the
#' \code{setDataPath} function.
#'
#'
#' @param type A character string, either \code{'static'} or \code{'yearly'}.
#' @param details A boolean. If \code{TRUE}, returns additional details about the available variables (sources, units, etc).
#'
#' @return A data frame.
#'
#' @examples
#'
#' \dontrun{
#' vars.df <- getVarNames(type = 'static', details = TRUE)
#' }
#'
#' @export
getVarNames <- function(type = c('static', 'yearly'), details = FALSE) {

  type <- match.arg(type)

  if (type == 'static') {
    staticPath <- getDataPath('static')
    staticHeader <- as.character(read.csv(file = staticPath, header = FALSE, nrows = 1, stringsAsFactors = FALSE)[1,])
    staticHeader.df <- data.frame(name = staticHeader)
    out <- staticHeader.df
  } else {
    yearlyPath <- getDataPath('yearly')
    yearlyHeader <- as.character(read.csv(file = yearlyPath, header = FALSE, nrows = 1, stringsAsFactors = FALSE)[1,])
    yearlyHeader.df <- data.frame(name = yearlyHeader)
    out <- yearlyHeader.df
  }

  if (details) {
    data(priogrid)
    out <- merge(x = out, y = priogrid.df, all.x = TRUE, all.y = FALSE, sort = FALSE)
    out <- out[order(out$id),]
  }

  return(out)
}

makeNumeric <- function(df) {
  for (i in 1:ncol(df)) {
    df[,i] <- as.numeric(df[,i])
  }
  return(df)
}

getDataFrame <- function(staticNames=c(), yearlyNames=c()) {

  out <- vector('list', 2)
  names(out) <- c('static', 'yearly')

  if (length(staticNames) > 0) {

    staticPath <- getDataPath('static')
    staticHeader.df <- getVarNames('static')
    selectVars <- c('gid', 'row', 'col', staticNames)
    selectVarsIdx <- which(staticHeader.df$name %in% selectVars)
    colClasses <- rep('NULL', nrow(staticHeader.df))
    colClasses[selectVarsIdx] <- NA
    static.df <- read.csv(staticPath, header = TRUE, colClasses = colClasses, stringsAsFactors = FALSE)
    static.df <- makeNumeric(static.df)
    out[['static']] <- static.df
  }

  if (length(yearlyNames) > 0) {

    yearlyPath <- getDataPath('yearly')
    yearlyHeader.df <- getVarNames('yearly')
    selectVars <- c('gid', 'year', yearlyNames)
    selectVarsIdx <- which(yearlyHeader.df$name %in% selectVars)
    colClasses <- rep('NULL', nrow(yearlyHeader.df))
    colClasses[selectVarsIdx] <- NA
    yearly.dt <- data.table::fread(yearlyPath, header = TRUE, colClasses = colClasses, stringsAsFactors = FALSE, showProgress = FALSE)
    yearly.df <- as.data.frame(yearly.dt)
    rm(yearly.dt)
    yearly.df <- makeNumeric(yearly.df)
    out[['yearly']] <- yearly.df
  }

  return(out)
}

getAvailYears <- function(year, value) {
  y <- unique(year[!is.na(value)])
  y <- y[order(y)]
  return(y)
}


#' Get PRIO GRID Data as a Raster
#'
#' Returns the requested PRIO GRID variables as a global raster at 0.5 x 0.5 decimal degree resolution.
#'
#' Note that only variables contained in the PRIO GRID csv files linked by the \code{setDataPath} function will
#' be available. Download the PRIO GRID csv files at \url{http://http://grid.prio.org/#/download}.
#'
#' @param names The names of the variables to be transformed into a raster. A warning will be issues for variables
#' that are unavailable.
#' @param years An integer vector indicating the years for which the yearly variables should be returned. Ignored for
#' static variables.
#'
#' @return A list containing (1) a VeloxRaster object and (2) a data frame with meta data about the PRIO GRID variables.
#'
#' @examples
#'
#' \dontrun{
#' prio.ls <- getPrioRaster(names = 'nlights_calib_mean', years = 2005:2010)
#' vx <- prio.ls$raster
#' meta.df <- prio.ls$meta
#' }
#'
#' @export
#' @useDynLib pgrid
#' @importFrom Rcpp evalCpp
getPrioRaster <- function(names, years=c()) {

  ## Get yearly variables

  # Empty containers (in case no yearly vars found)
  yearlyMatrix.ls <- list()
  yearlyRaster.df <- data.frame(name = c(), year = c(), stringsAsFactors = FALSE)

  if (isDataPathAvail('yearly')) {
    yearlyVars.df <- getVarNames(type = 'yearly')
    if (any(names %in% yearlyVars.df$name)) {

      # Prepare vectors for storing raster meta info
      yearlyRasterNames <- c()
      yearlyRasterYears <- c()

      # Get raster data as DF
      yearlyNames <- names[names %in% yearlyVars.df$name]
      yearlyData.df <- getDataFrame(yearlyNames = yearlyNames)$yearly

      # Iterate over variables
      yearlyMatrix.ls <- vector('list', length(yearlyNames))
      for (i in 1:length(yearlyNames)) {

        thisName <- yearlyNames[i]

        # Get relevant data vectors
        thisCol <- which(names(yearlyData.df) == thisName)
        yearVec <- yearlyData.df$year
        valVec <- yearlyData.df[,thisCol]
        gidVec <- yearlyData.df$gid

        # Get relevant years
        yearsAvail <- getAvailYears(yearVec, valVec)
        if (length(years) > 0) {
          yearsSelect <- yearsAvail[yearsAvail %in% years]
        } else {
          yearsSelect <- yearsAvail
        }

        # Skip to next variable if no relevant years
        if (length(yearsSelect) == 0) {
          next
        }

        # Get matrix for each year
        variableMatrix.ls <- vector('list', length(yearsSelect))
        for (j in 1:length(yearsSelect)) {

          thisYear <- yearsSelect[j]

          thisYearVal <- valVec[yearVec == thisYear]
          thisYearGid <- gidVec[yearVec == thisYear]
          thisMat <- getPrioMatrix(thisYearGid, thisYearVal)
          variableMatrix.ls[[j]] <- thisMat

          yearlyRasterNames <- c(yearlyRasterNames, thisName)
          yearlyRasterYears <- c(yearlyRasterYears, thisYear)
        }
        yearlyMatrix.ls[[i]] <- variableMatrix.ls
      }

      yearlyMatrix.ls <- unlist(yearlyMatrix.ls, recursive = FALSE)  # Flatten list
      yearlyRaster.df <- data.frame(name = yearlyRasterNames, year = yearlyRasterYears, stringsAsFactors = FALSE)
      names <- names[!(names %in% yearlyNames)]
    }
  }


  ## Get static variables

  # Empty containers (in case no static vars)
  staticMatrix.ls <- list()
  staticRaster.df <- data.frame(name = c(), year = c(), stringsAsFactors = FALSE)

  if (isDataPathAvail('static')) {
    staticVars.df <- getVarNames(type = 'static')
    if (any(names %in% staticVars.df$name)) {

      # Prepare vectors for storing raster meta info
      staticRasterNames <- c()
      staticRasterYears <- c()

      # Get raster data as DF
      staticNames <- names[names %in% staticVars.df$name]
      staticData.df <- getDataFrame(staticNames = staticNames)$static

      # Iterate over variables
      staticMatrix.ls <- vector('list', length(staticNames))
      for (i in 1:length(staticNames)) {

        thisName <- staticNames[i]

        # Get relevant data vectors
        thisCol <- which(names(staticData.df) == thisName)
        valVec <- staticData.df[,thisCol]
        gidVec <- staticData.df$gid

        # Get matrix
        thisMat <- getPrioMatrix(gidVec, valVec)
        staticMatrix.ls[[i]] <- thisMat

        # Store meta info
        staticRasterNames <- c(staticRasterNames, thisName)
        staticRasterYears <- c(staticRasterYears, NA)
      }

      staticRaster.df <- data.frame(name = staticRasterNames, year = staticRasterYears, stringsAsFactors = FALSE)
      names <- names[!(names %in% staticNames)]
    }
  }

  ## Issue warning for unfound vars
  orphanVars <- names
  if (length(orphanVars) > 0) {
    warning(paste('Variable(s)', paste(orphanVars, collapse = ','), 'not found in linked PRIO GRID csv files.'))
  }

  ## Combine meta data
  meta.df <- rbind(yearlyRaster.df, staticRaster.df)

  ## Exit early if no raster data found
  if (nrow(meta.df) == 0) {
    out.ls <- list(raster = NULL, meta = meta.df)
    return(out.ls)
  }

  ## Add band info to meta data
  meta.df$band <- 1:nrow(meta.df)

  ## Combine matrices, make VeloxRaster
  vx.ls <- c(yearlyMatrix.ls, staticMatrix.ls)
  vx.nms <- c(paste(yearlyRaster.df$name, yearlyRaster.df$year, sep = '_'), staticRaster.df$name)
  names(vx.ls) <- vx.nms
  vx <- velox::velox(x = vx.ls, extent = c(-180, 180, -90, 90),
              res = c(0.5, 0.5), crs = '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ')

  out.ls <- list(raster = vx, meta = meta.df)
  return(out.ls)
}





