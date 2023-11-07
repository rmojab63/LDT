## code to prepare `data.wdi` dataset

# Specify the location where the raw data is saved:
RawDataFolder <- "D:/Data/WDI"

#' Loads and Aggregates WDI Data
#'
#' Use this function to aggregate and reshape data from the World Development Indicators dataset. It creates a 'Country-Series' table.
#'
#' @param dirPath A character string representing the path to the WDI dataset directory. It should contain the following files: ‘WDICountry-Series.csv’, ‘WDIData.csv’, ‘WDICountry.csv’, and ‘WDISeries.csv’. These files can be downloaded from the WDI website.
#' @param minYear An integer representing the starting year for data aggregation.
#' @param maxYear An integer representing the ending year for data aggregation.
#' @param aggregateFun A function for aggregation with the following arguments: \code{data}, \code{code}, \code{name}, \code{unit}, \code{definition}, \code{aggMethod}, where \code{data} is the data-points from \code{minYear} to \code{maxYear}, \code{unit} is the unit of measurement, \code{definition} is the long definition of the series, and \code{aggMethod} is the method of aggregation.
#' @param keepFun A function that determines whether to keep or omit columns of the final data matrix. It has a \code{data} argument which is the data of the column. It can check the variance or count the number of available data points and omit a variable from the analysis.
#'
#' @return A list with the following items:
#' \item{result}{A matrix with countries as rows and series as columns. Each data point represents the aggregated value from \code{minYear} to \code{maxYear}.}
#' \item{countries}{A list containing information about the countries.}
#' \item{series}{A list containing information about various series.}
#'
#' @export
#'
#' @examples
#' # This example requires external data. Download the data and run it:
#'
#' \donttest{
#' try({
#'    path.dir <- "D:/Data/WDI" # This path must be valid
#'    # Define a function that calculates the long-run growth rates:
#'    aggregateFun <- function(data, code, name, unit, definition, aggMethod) {
#'      isPerc <- unit == "%" || grepl(".ZG", code)
#'      if (isPerc) NA
#'      else get.longrun.growth(data, FALSE, TRUE, 30, 5, isPerc)
#'    }
#'    # Add some rules for removing variables from the analysis:
#'    keepFun <- function(data) {
#'         var(data, na.rm = TRUE) > 1e-12 && sum((is.na(data)) == FALSE) >= 50
#'    }
#'    data <- data.wdi.agg(path.dir, 1960, 2020, aggregateFun, keepFun)
#' })
#' }
#'
#' @importFrom utils read.csv
#' @importFrom stats var
#'
data.wdi.agg <- function(dirPath, minYear = 1960, maxYear = 2020, aggregateFun = NULL, keepFun = NULL) {
  cal <- match.call()
  # read countries
  countryFileName <- "WDICountry.csv"
  raw_c <- read.csv(paste0(dirPath, "/", countryFileName))
  countriesCodes <- raw_c$Country.Code

  # read series
  seriesFileName <- "WDISeries.csv"
  raw_s <- read.csv(paste0(dirPath, "/", seriesFileName))
  seriesCodes <- raw_s$Series.Code

  # final data matrix
  result <- matrix(NA, length(countriesCodes), length(seriesCodes),
                   dimnames = list(countriesCodes, seriesCodes)
  )

  # meta: WDIData.csv (aggregate data)
  dataFileName <- "WDIData.csv"
  indCountry <- 1
  indCountryCode <- 2
  indSeries <- 3
  indSeriesCode <- 4
  indData <- 5
  startYear <- 1960

  firstInd <- minYear - startYear + indData
  lastInd <- maxYear - startYear + indData

  if (firstInd < indData || lastInd < indData) {
    stop(paste0("invalid date. 'data' starts at ", startYear))
  }

  #print("reading CSV file. This might be time-consuming ...")
  con <- file(paste0(dirPath, "/", dataFileName), open = "rt")

  # all <- length(readLines(con)) I'm not sure about the efficiency of
  # this and therefore, I do not add a counter

  ind <- 0
  while (TRUE) {
    ind <- ind + 1
    cells <- strsplit(readLines(con, n = 1), "\",\"")
    if (length(cells) == 0) {
      break()
    }
    if (ind == 1) {
      next()
    }
    cells <- cells[[1]]
    i <- which(countriesCodes == cells[indCountryCode])
    if (length(i) > 0) {
      j_code <- cells[[indSeriesCode]]
      j <- which(seriesCodes == j_code)
      j_unit <- raw_s$Unit.of.measure[[j]]
      j_name <- raw_s$Indicator.Name[[j]]
      j_definition <- raw_s$Long.definition[[j]]
      j_aggMethod <- raw_s$Aggregation.method
      if (length(j) > 0) {
        result[[i[[1]], j[[1]]]] <- aggregateFun(
          as.numeric(cells[firstInd:lastInd]), j_code,
          j_name, j_unit, j_definition, j_aggMethod
        )
      }
    }
  }
  close(con)
  #print("reading CSV file finished.")

  # NA columns
  c_inds <- which(as.logical(sapply(as.data.frame(result), function(data) keepFun(data))))
  result <- result[, c_inds]

  # don't remove any rows (for example empty rows) or you
  # have to adjust the data from different calls

  list(
    data = result,
    countries = list(
      code = countriesCodes,
      name = raw_c$Short.Name,
      isCountry = raw_c$System.of.National.Accounts != "" # assuming that 'System of National
      #                                                     Accounts' is not empty for countries
    ),
    series = list(
      codes = seriesCodes[c_inds],
      names = raw_s$Indicator.Name[c_inds],
      descriptions = raw_s$Long.definition[c_inds],
      units = raw_s$Unit.of.measure[c_inds],
      topics = raw_s$Topic[c_inds],
      aggregationMethods = raw_s$Aggregation.method[c_inds]
    ),
    call = cal
  )
}

has_all <- function(keywords, text) {
  if (length(text) == 0)
    return(FALSE)
  if (length(keywords) == 0)
    stop(paste0("invalid keywords. It is empty. text='", text,"'"))
  for (key in keywords){
    if (length(key) == 0)
      next()
    if (grepl(key, text, fixed = TRUE) == FALSE)
      return(FALSE)
  }
  return(TRUE)
}


#' Finds a Series in WDI Dataset
#'
#' Use this function to search for a series in the WDI dataset by its name, description, and some other attributes.
#'
#' @param series This must be the \code{series} field in the output of the \code{\link{data.wdi.agg}} function.
#' @param keywords A character array representing the keywords to be used for searching.
#' @param searchName A logical value indicating whether to search in the names field or not. If \code{FALSE}, it does not search in the names field.
#' @param searchDesc A logical value indicating whether to search in the description or not. If \code{FALSE}, it does not search in the description.
#' @param topicKeywords A character array representing the keywords that must be contained in the topic of a matched case.
#' @param findOne A logical value indicating whether to raise an error if more than one series is found. The default is \code{FALSE}.
#'
#' @return If \code{findOne} is \code{TRUE}, it returns a series information. Otherwise, it returns a list with series information.
#' @export
#'
#' @examples
#' #data <- data.wdi.agg() # this is time-consuming and requires WDI dataset files
#' #res <- data.wdi.search(data$series, c("GDP per capita"),
#' #                        TRUE, topicKeywords = "national account")
data.wdi.search <- function(series, keywords, searchName = TRUE,
                            searchDesc = FALSE, topicKeywords = NULL, findOne = FALSE) {
  if (is.null(series)) {
    stop("series is null. use Data_Wdi() function to generate it.")
  }
  if (is.vector(keywords) == FALSE) {
    keywords <- c(keywords)
  }
  keywords <- tolower(keywords)
  topicKeywords <- tolower(topicKeywords)

  res <- list()
  for (i in c(1:length(series$codes))) {
    code <- tolower(series$codes[[i]])
    name <- tolower(series$names[[i]])
    description <- tolower(series$descriptions[[i]])
    topic <- tolower(series$topics[[i]])

    if (keywords[[1]] == code ||
        (searchName && has_all(keywords, name)) ||
        (searchDesc && has_all(keywords, description))) {
      if (length(topicKeywords) == 0 || has_all(topicKeywords, topic)) {
        res[[length(res) + 1]] <- list(
          code = series$codes[[i]], # don't use lowercase
          name = series$names[[i]],
          description = series$descriptions[[i]],
          unit = series$units[[i]],
          topic = series$topics[[i]],
          aggregationMethod = series$aggregationMethods[[i]]
        )
      }
    }
  }
  if (length(res) == 0) {
    stop("WDI series not found.")
  }
  if (findOne) {
    if (length(res) > 1) {
      stop("More than one WDI series is found.")
    }
    return(res[[1]])
  } else {
    return(res)
  }
}



# We are going to study the effect of variables before the 'splitYear'  GDP per capita growth after this year
splitYear <- 2005

# Aggregate data before 'splitYear', using 'tdata' package:
wdi_lr_2005 <- tdata::data.wdi.agg(
  dirPath = RawDataFolder,
  minYear = 1960, maxYear = splitYear,
  aggregateFun = function(data, code, name, unit, definition, aggMethod) {
    isPerc <- unit == "%" || grepl(".ZG", code)
    if (isPerc) NA else tdata::get.longrun.growth(data, TRUE, isPerc, 20, 2, TRUE)
  },
  keepFun = function(X) {
    var(X, na.rm = TRUE) > 1e-12 && sum((is.na(X)) == FALSE) >= 50
  }
)

# Aggregate data after 'splitYear', using 'tdata' package:
wdi_lr_2006_ <- tdata::data.wdi.agg(
  dirPath = RawDataFolder,
  minYear = splitYear+1, maxYear = 2021,
  aggregateFun = function(data, code, name, unit, definition, aggMethod) {
    isPerc <- unit == "%" || grepl(".ZG", code)
    if (isPerc) NA else tdata::get.longrun.growth(data, TRUE, isPerc, 2, 4, TRUE)
  },
  keepFun = function(X) {
    var(X, na.rm = TRUE) > 1e-12 && sum((is.na(X)) == FALSE) >= 50
  }
)

# Create the data matrix with variables in the columns
#     by combining the dependent after splitYear and regressors before that year.
serCode1 <- tdata::data.wdi.search(wdi_lr_2005$series, c("GDP", "per capita", "constant", "us"))[[1]]$code
#serCode2 <- tdata::data.wdi.search(wdi_lr_2005$series, c("consumer price index"))[[1]]$code
y <- wdi_lr_2006_$data[, which(colnames(wdi_lr_2006_$data) %in% c(serCode1#, serCode2
)), drop = FALSE]
x <- as.data.frame(wdi_lr_2005$data)

# Remove highly correlated columns and non-country observations:
x <- x[, which(sapply(x, function(v) abs(cor(y[,1], v, use = "complete.obs")) < 0.999)), drop = FALSE]
#x <- x[, which(sapply(x, function(v) abs(cor(y[,2], v, use = "complete.obs")) < 0.999)), drop = FALSE]
wdi_y <- y[wdi_lr_2005$countries$isCountry, , drop = FALSE]
wdi_y[is.nan(wdi_y)] <- NA

x <- x[wdi_lr_2005$countries$isCountry, , drop = FALSE]
if (any(rownames(x) != names(y))) {
  stop("Invalid operation.")
}

# rearrange data such that the lag of the dependent variable
#    is in the first column
ind <- which(colnames(x) %in% c(serCode1))
cn <- colnames(x)
x <- cbind(x[, ind, drop = F], x[, 1:(ind - 1), drop = F], x[, ((ind + 1):ncol(x)), drop = F])
colnames(x) <- c(cn[ind], cn[1:(ind-1)], cn[(ind+1):length(cn)])

# add intercept:
wdi_x <- as.matrix(x)

# make sure they are number:
wdi_x <- matrix(as.numeric(wdi_x), ncol = ncol(wdi_x), dimnames = list(NULL, colnames(wdi_x)))
wdi_x[is.nan(wdi_x)] <- NA

# get names
x_names <- lapply(colnames(wdi_x), function(c){
  ind <- which(wdi_lr_2005$series$codes == c)
  if (length(ind) > 0)
    return(list(code = c, name = wdi_lr_2005$series$names[[ind[[1]]]]))
  return(list(code = c, name = NULL))
})

# A this point, `wdi_growth` is a table with countries in the rows and variables in the columns.
data.wdi <- list(y = wdi_y, x = wdi_x, splitYear = splitYear, names = x_names)



usethis::use_data(data.wdi, overwrite = TRUE)
