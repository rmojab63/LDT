

#' Loads and Aggregates WDI Data
#'
#' Use it to aggregate and reshape data from the World Development Indicators dataset. It creates Creates a 'Country-Series' table.
#'
#' @param dirPath (character) Path to the WDI dataset directory. It should contain the following files: ‘WDICountry-Series.csv’, ‘WDIData.csv’, ‘WDICountry.csv’, and ‘WDISeries.csv’. These files can be downloaded from the WDI website.
#' @param minYear (integer) Starting year for data aggregation.
#' @param maxYear (integer) Ending year for data aggregation.
#' @param aggregateFun (function) Aggregation function with the following arguments: \code{data}, \code{code}, \code{name}, \code{unit}, \code{definition}, \code{aggMethod}, where \code{data} is the data-points from \code{minYear} to \code{maxYear}, \code{unit} is the unit of measurement, \code{definition} is the long definition of the series, \code{aggMethod} is the method of aggregation.
#' @param keepFun (function) A Function that determines whether to keep or omit columns of the final data matrix. It has \code{data} argument which is the data of the column. It can check the variance or count the number of available data points and omit a variable from the analysis.
#'
#' @return A list with the following arguments:
#' \tabular{ll}{
#' \code{result} \tab A matrix with countries as rows and series as columns. Each data point represents the aggregated value from \code{minYear} to \code{maxYear}. \cr
#' \code{countries} \tab A list containing information about the countries. \cr
#' \code{series} \tab A list containing information about various series.
#' }
#' @export
#'
#' @examples
#' # This example requires external data. Download the data and run them:
#'
#' \donttest{
#' path.dir <- "D:/Data/WDI" # This path must be valid
#' # Define a function that calculates the long-run growth rates:
#' aggregateFun <- function(data, code, name, unit, definition, aggMethod) {
#'   isPerc <- unit == "%" || grepl(".ZG", code)
#'   if (isPerc) NA else LongrunGrowth(data, 30, 5, FALSE, TRUE, isPerc)
#' }
#' # Add some rules for removing the variable from the analysis:
#' keepFun <- function(data) {
#'      var(data, na.rm = TRUE) > 1e-12 && sum((is.na(data)) == FALSE) >= 50
#' }
#' data <- data.wdi.agg(path.dir, 1960, 2020, aggregateFun, keepFun)
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
    stop(paste0("Invalid date. data starts at ", startYear))
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
    stop(paste0("Invalid keywords. It is empty. text='", text,"'"))
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
#' Use it to search for a series in the WDI dataset by its name, description, and some other attributes.
#'
#' @param series This must be the \code{series} field in the output of the \code{\link{data.wdi.agg}} function.
#' @param keywords (character array) Keywords to be used for searching.
#' @param searchName If \code{FALSE}, it does not search in the names field.
#' @param searchDesc If \code{FALSE}, it does not search in the description
#' @param topicKeywords If given, topic of a matched case must contain these keywords, too.
#' @param findOne Raises error if \code{TRUE} and more than 1 series is found. default is \code{FALSE}.
#'
#' @return If \code{findOne} is \code{TRUE}, it is a series information. Otherwise, it is a list with series information.
#' @export
#'
#' @examples
#' #data <- data.wdi.agg() # this is time-consuming and requires WDI dataset files
#' #res <- data.wdi.search(data$series, c("GDP per capita"),
#' #                        TRUE, topicKeywords = "national account")
#'
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


#' private method to generate dummy variable
#'
#' @param table data
#' @param colName categorical column
#' @param pre a string to put before the name of the variables
#' @param min_unique_skip if number of unique values is equal or larger, it returns \code{NULL}
#' @param uniques if not \code{NULL}, it skips finding unique values and uses
#' the given list. Also, if \code{colName} column is missing,
#' it creates zero variables for the given items
#'
#' @return \code{data} (list of dummy variables) and \code{uniques} (unique values)
getDummy <- function(table, colName, pre = "",
                     min_unique_skip = 6, uniques = NULL) {
  col <- tryCatch(table[, colName], error = function(e) NULL)
  if (is.null(col) && is.null(uniques)) {
    stop(paste0("column is missing. colName = ", colName))
  }

  if (is.null(uniques)) {
    uniques <- c(unique(col))
  }
  if (length(uniques) >= min_unique_skip) {
    return(NULL)
  }

  ns <- paste0(pre, colName, "_", sapply(uniques, function(u) gsub(" ", "_", u)))

  res <- lapply(uniques, function(v) {
    if (is.null(col)) {
      numeric(nrow(table))
    } else {
      as.numeric(col == v)
    }
  })
  names(res) <- ns

  list(data = res, uniques = uniques)
}


#' Load 'Berka' Dataset
#'
#' Use it to load and combine tables from the ‘Berka’ dataset to create a unified data table.
#'
#' @param dirPath Path to the downloaded data directory.
#' @param positive Determines the positive class. There are four types of loans:
#' 'A' stands for contract finished, no problems, 'B' stands for contract finished,
#' loan not payed, 'C' stands for running contract, OK so far, 'D' stands
#' for running contract, client in debt.
#' @param negative Similar to \code{positive} for negative class.
#' @param rateFun A function to calculate interest rate in loans. Its arguments are: \code{amount}, \code{duration}, \code{paymentPerMonth}.
#'
#' @return A data.frame with the following columns:
#' \tabular{ll}{
#' \code{loan_id} \tab record identifier \cr
#' \code{status} \tab original status of the data (A, B, C, or D) \cr
#' \code{label} \tab status of paying off the loan transformed to numeric (0,1) by using \code{positive} and \code{negative} arguments. value=1 means default. \cr
#' \code{amount} \tab amount of money \cr
#' \code{payments} \tab monthly payments \cr
#' \code{rate} \tab rates calculated by \code{rateFun} function \cr
#' \code{duration_# (#=12,24,36,48,60)} \tab dummy variables for the duration of the loan \cr
#' \code{account_frequency_?} \tab dummy variables for the frequency of issuance of statements. ?="POPLATEK MESICNE" stands for monthly issuance, ?="POPLATEK TYDNE" stands for weekly issuance, ?="POPLATEK PO OBRATU" stands for issuance after transaction \cr
#' \code{order_num} \tab number of the payment orders issued for the account of the loan \cr
#' \code{order_sum_amount} \tab sum of amounts of the payment orders issued for the account of the loan \cr
#' \code{order_related_account_num} \tab unique number of 'account of the recipient' in the payment orders issued for the account of the loan \cr
#' \code{order_related_bank_num} \tab unique number of 'bank of the recipient' in the payment orders issued for the account of the loan \cr
#' \code{order_has_?} \tab dummy variables fo 'characterization of the payment' in the payment orders issued for the account of the loan \cr
#' \code{trans_?num} \tab number of transactions dealt with the account of the loan (in different groups) \cr
#' \code{trans_?amount_mean} \tab mean of 'amount of money' in the transactions dealt with the account of the loan (in different groups) \cr
#' \code{trans_?amount_div_balance} \tab mean of 'amount of money'/'balance after transaction' in the transactions dealt with the account of the loan (in different groups) \cr
#' \code{trans_related_account_num} \tab unique number of 'account of the partner' in the transactions dealt with the account of the loan \cr
#' \code{trans_related_account_num} \tab unique number of 'bank of the partner' in the transactions dealt with the account of the loan \cr
#' \code{dist_inhabitants_num} \tab no. of inhabitants in the location of the branch of the account of the loan \cr
#' \code{dist_muni_#1#2} \tab no. of municipalities with inhabitants #1-#2 in the location of the branch of the account of the loan \cr
#' \code{dist_cities_num} \tab no. of cities in the location of the branch of the account of the loan \cr
#' \code{dist_ratio_urban_inhabitants} \tab ratio of urban inhabitants in the location of the branch of the account of the loan \cr
#' \code{dist_avg_salary} \tab average salary in the location of the branch of the account of the loan \cr
#' \code{dist_unemployment95} \tab unemployment rate '95 in the location of the branch of the account of the loan \cr
#' \code{dist_unemployment96} \tab unemployment rate '96 in the location of the branch of the account of the loan \cr
#' \code{dist_entrepreneurs_num_per1000} \tab no. of entrepreneurs per 1000 inhabitants in the location of the branch of the account of the loan \cr
#' \code{dist_crimes95_num} \tab no. of committed crimes '95 in the location of the branch of the account of the loan \cr
#' \code{dist_crimes96_num} \tab no. of committed crimes '96 in the location of the branch of the account of the loan
#' }
#'
#' @export
#' @importFrom utils read.csv
data.berka.loan <- function(dirPath,
                           positive = c("B", "D"), negative = c("A", "C"),
                           rateFun = function(amount, duration, paymentPerMonth) {
                             ((paymentPerMonth * duration) / amount - 1) * 100
                           }) {
  loans0 <- as.data.frame(read.csv(paste0(dirPath, "/loan.asc"), header = TRUE, sep = ";"))
  # get labels
  loans <- data.frame(loan_id = loans0$loan_id)
  loans$status <- loans0$status
  loans$label <- as.numeric(sapply(
    c(loans0$status),
    function(s) {
      if (any(positive == s)) 1 else if (any(negative == s)) 0 else NA
    }
  ))

  loans$amount <- loans0$amount
  loans$duration <- loans0$duration
  loans$payments <- loans0$payments
  loans$account_id <- loans0$account_id

  inv_na <- which(is.na(loans$label))
  if (length(inv_na) > 0) {
    loans <- loans[-inv_na, ]
  } # remove invalid rows

  account_id <- loans$account_id
  loans$account_id <- NULL

  loans <- cbind(loans, data.frame(rate = rateFun(loans$amount, loans$duration, loans$payments)))
  loans <- cbind(
    loans[, -c(which(colnames(loans) == "duration"))],
    getDummy(loans, "duration")$data
  ) # convert to dummy


  # ACCOUNT
  account <- as.data.frame(read.csv(paste0(dirPath, "/account.asc"), header = TRUE, sep = ";"))
  account <- account[c(sapply(
    account_id,
    function(id) which(account$account_id == id)
  )), ] # they are unique
  loans <- cbind(
    loans,
    getDummy(account, "frequency", "account_")$data
  ) # frequency of issuance of statements:
  #                                                   Monthly, Weekly or issuance after transaction
  district_id <- account$district_id


  # ORDER
  order <- read.csv(paste0(dirPath, "/order.asc"), header = TRUE, sep = ";")
  order_inds <- sapply(account_id, function(id) which(order$account_id == id))
  loans$order_num <- as.numeric(sapply(
    order_inds,
    function(inds) length(inds)
  )) # number of orders with the account
  loans$order_sum_amount <- as.numeric(sapply(
    order_inds,
    function(inds) sum(order$amount[inds])
  )) # sum of ordered amounts
  loans$order_related_account_num <- as.numeric(sapply(
    order_inds,
    function(inds) length(unique(order$account_to[inds]))
  )) # number of related accounts
  loans$order_related_bank_num <- as.numeric(sapply(
    order_inds,
    function(inds) length(unique(order$bank_to[inds]))
  )) # number of related accounts
  loans$order_has_insurrance <- as.numeric(sapply(
    order_inds,
    function(inds) any(order$k_symbol[inds] == "POJISTNE")
  ))
  loans$order_has_household <- as.numeric(sapply(
    order_inds,
    function(inds) any(order$k_symbol[inds] == "SIPO")
  ))
  loans$order_has_leasing <- as.numeric(sapply(
    order_inds,
    function(inds) any(order$k_symbol[inds] == "LEASING")
  ))
  loans$order_has_loan <- as.numeric(sapply(
    order_inds,
    function(inds) any(order$k_symbol[inds] == "UVER")
  ))

  # TRANSACTIONS
  trans <- read.csv(paste0(dirPath, "/trans.asc"), header = TRUE, sep = ";")
  trans <- trans[c(sapply(
    c(1:nrow(trans)),
    function(i) any(account_id == trans$account_id[[i]])
  )), ] # remove irrelevant
  trans_inds <- sapply(account_id, function(id) which(trans$account_id == id))
  prenames <- c(
    "", "type_credit_", "type_withdrawal_",
    "op_credit_withdrawal_", "op_credit_cash_", "op_collect_bank_",
    "op_cash_withdrawal_", "op_remitence_bank_"
  )
  i <- 0
  for (c in list(
    trans_inds, # filter this for more variables...

    sapply(trans_inds, function(inds) which(trans$type[inds] == "PRIJEM")),
    sapply(trans_inds, function(inds) which(trans$type[inds] == "VYDAJ")),
    sapply(trans_inds, function(inds) which(trans$operation[inds] == "VYBER KARTOU")),
    sapply(trans_inds, function(inds) which(trans$operation[inds] == "VKLAD")),
    sapply(trans_inds, function(inds) which(trans$operation[inds] == "PREVOD Z UCTU")),
    sapply(trans_inds, function(inds) which(trans$operation[inds] == "VYBER")),
    sapply(trans_inds, function(inds) which(trans$operation[inds] == "PREVOD NA UCET"))
  )) {
    i <- i + 1
    loans[, paste0("trans_", prenames[[i]], "num")] <- as.numeric(sapply(
      c,
      function(inds) length(inds)
    ))
    loans[, paste0("trans_", prenames[[i]], "amount_mean")] <- as.numeric(sapply(
      c,
      function(inds) mean(trans$amount[inds])
    ))
    loans[, paste0("trans_", prenames[[i]], "amount_div_balance")] <- as.numeric(sapply(
      c,
      function(inds) mean(trans$amount[inds] / trans$balance[inds])
    ))
  }
  loans$trans_related_account_num <- as.numeric(sapply(
    trans_inds,
    function(inds) length(unique(trans$bank_to[inds]))
  )) # number of related accounts
  loans$trans_related_bank_num <- as.numeric(sapply(
    trans_inds,
    function(inds) length(unique(order$account_to[inds]))
  )) # number of related accounts


  # DEMOGRAPHIC
  district <- read.csv(paste0(dirPath, "/district.asc"), header = TRUE, sep = ";")
  district_inds <- sapply(district_id, function(id) which(district$A1 == id))
  loans$dist_inhabitants_num <- as.numeric(sapply(
    district_inds,
    function(inds) district$A4[[inds]]
  ))
  loans$dist_muni_0_499 <- as.numeric(sapply(
    district_inds,
    function(inds) district$A5[[inds]]
  )) # no. of municipalities with inhabitants < 499
  loans$dist_muni_500_1999 <- as.numeric(sapply(
    district_inds,
    function(inds) district$A6[[inds]]
  )) # no. of municipalities with inhabitants 500-1999
  loans$dist_muni_2000_9999 <- as.numeric(sapply(
    district_inds,
    function(inds) district$A7[[inds]]
  )) # no. of municipalities with inhabitants 2000-9999
  loans$dist_muni_100000_ <- as.numeric(sapply(
    district_inds,
    function(inds) district$A8[[inds]]
  )) # no. of municipalities with inhabitants >10000
  loans$dist_cities_num <- as.numeric(sapply(
    district_inds,
    function(inds) district$A9[[inds]]
  ))
  loans$dist_ratio_urban_inhabitants <- as.numeric(sapply(
    district_inds,
    function(inds) district$A10[[inds]]
  ))
  loans$dist_avg_salary <- as.numeric(sapply(
    district_inds,
    function(inds) district$A11[[inds]]
  ))
  loans$dist_unemployment95 <- as.numeric(sapply(
    district_inds,
    function(inds) district$A12[[inds]]
  ))
  loans$dist_unemployment96 <- as.numeric(sapply(
    district_inds,
    function(inds) district$A13[[inds]]
  ))
  loans$dist_entrepreneurs_num_per1000 <- as.numeric(sapply(
    district_inds,
    function(inds) district$A14[[inds]]
  ))
  loans$dist_crimes95_num <- as.numeric(sapply(
    district_inds,
    function(inds) district$A15[[inds]]
  ))
  loans$dist_crimes96_num <- as.numeric(sapply(
    district_inds,
    function(inds) district$A16[[inds]]
  ))

  # Client
  # Unfortunately there is no account_id in this table. just a distrinct_id
  # client = read.csv(paste0(dirPath,"/client.asc"), header = TRUE, sep = ";")
  # client_inds = sapply(district_id, function(id)which(client$district_id == id))
  # loans$client_birth_date =  ...

  # disposition
  # this is irrelevant : "only owner can issue permanent orders and ask for a loan"
  # disps = read.csv(paste0(dirPath,"/disp.asc"), header = TRUE, sep = ";")
  # disps_inds = sapply(account_id, function(id)which(disps$account_id == id))
  # unique(lengths(disps_inds)) !!! some has two indexes (it seems that all are owner w)
  # loans$disp_isowner = as.numeric(sapply(disps_inds,
  #     function(inds)any(disps$type[inds]=="OWNER")))

  # CARDS
  # cards = read.csv(paste0(dirPath,"/card.asc"), header = TRUE, sep = ";")
  # cards_inds = sapply(account_id, function(id)which(cards$account_id == id))
  # loans$card_num = as.numeric(sapply(cards_inds, function(inds)length(inds)))


  loans
}



#' Load 'Vesta' Dataset
#'
#' Use it to create fraud-series table from 'Vesta' dataset (aka IEEE-CIS Fraud Detection).
#'
#' @param dirPath Path to the downloaded data directory.
#' @param training If \code{FALSE}, it loads test data.
#' @param tranDumCols A list with \code{name} and \code{values} of (categorical) columns in 'transaction' file to be converted to dummy variables. If \code{training} is \code{FALSE} and this is \code{NULL}, a warning is raised.
#' @param idenDumCols Similar to \code{tranDumCols} but for 'identity' file.
#' @param catMinSkip If \code{tranDumCols} or \code{idenDumCols} is \code{NULL}, for a categorical variable, if number of unique values is equal or larger than this value, it is omitted.
#'
#' @return A list with the following items:
#' \tabular{ll}{
#' \code{data} \tab A \code{data.frame} with the data. \cr
#' \code{tranDumCols} \tab A list with \code{name} and \code{values} in 'transaction' data, used for creating the dummy variable. \cr
#' \code{idenDumCols} \tab A list with \code{name} and \code{values} in 'identity' data, used for creating the dummy variable.
#' }
#'
#' @export
#' @importFrom utils read.csv
data.vesta.fraud <- function(dirPath, training = TRUE,
                            tranDumCols = NULL, idenDumCols = NULL,
                            catMinSkip = 6) {
  # TODO: add an argument to group some items in categorical data
  #      for example, one might want to create the following groups for
  # P_emaildomain: 'gmail', 'yahoo', 'others'

  if (training == FALSE && (is.null(tranDumCols) || is.null(idenDumCols))) {
    warning("'?_dumCols' is not provided. The result might not
    be consistent with the training data.")
  }

  pre <- if (training) "train_" else "test_"

  transaction <- as.data.frame(
    read.csv(paste0(dirPath, "/", pre, "transaction.csv"),
             header = TRUE, sep = ","
    )
  )

  givenNames <- FALSE
  if (is.null(tranDumCols)) {
    tranDumCols <- c(
      "ProductCD", paste0("card", c(1:6)), "addr1", "addr2",
      "P_emaildomain", "R_emaildomain",
      paste0("M", c(1:9))
    )
  } else {
    givenNames <- TRUE
  }

  t_d <- list()
  tranDumCols0 <- list()
  for (tc in tranDumCols) {
    name <- if (givenNames) tc$name else tc
    uniques <- if (givenNames) tc$values else NULL
    res <- getDummy(transaction,
                    colName = name, "", catMinSkip,
                    uniques = uniques
    ) # note that in 'test', some columns might be missing. They are actually all 0
    if (is.null(res) == FALSE) {
      tranDumCols0[[length(tranDumCols0) + 1]] <- list(name = name, values = res$uniques)
      t_d <- append(t_d, res$data)
    }
    transaction[name] <- NULL
  }
  transaction <- cbind(transaction, t_d)


  # append identity
  identity <- as.data.frame(read.csv(paste0(dirPath, "/", pre, "identity.csv"),
                                     header = TRUE, sep = ","
  ))

  if (is.null(idenDumCols)) {
    idenDumCols <- c(
      "DeviceType", "DeviceInfo",
      paste0("id_", c(12:38))
    )
  }

  i_d <- list()
  idenDumCols0 <- list()
  for (tc in idenDumCols) {
    name <- if (givenNames) tc$name else tc
    uniques <- if (givenNames) tc$values else NULL
    res <- getDummy(identity,
                    colName = name, "", catMinSkip,
                    uniques = uniques
    ) # note that in 'test', some columns might be missing. They are actually all 0
    if (is.null(res) == FALSE) {
      idenDumCols0[[length(idenDumCols0) + 1]] <- list(name = name, values = res$uniques)
      i_d <- append(i_d, res$data)
    }
    identity[name] <- NULL
  }
  identity <- cbind(identity, i_d)

  inds <- sapply(transaction$TransactionID, function(id) match(id, identity$TransactionID))

  transaction <- cbind(transaction, identity[inds, ])

  return(list(data = transaction, tranDumCols = tranDumCols0, idenDumCols = idenDumCols0))
}


#' Load PCP Dataset
#'
#' Use it to load the 'IMF's Primary Commodity Prices' dataset and create Date-Series table.
#'
#' @param dirPath path to the downloaded data data.
#' @param makeReal uses the first column (which must be US-CPI) and converts nominal variables to real. If \code{TRUE}, \code{dirPath} must also contain a file with the US CPI.
#'
#' @return A list with data, descriptions, etc.
#'
#' @export
data.pcp <- function(dirPath, makeReal = FALSE) {
  fs <- list.files(dirPath, "*.xls")
  data0 <- readxl::read_excel(file.path(dirPath, fs[[1]]))

  descriptions <- c(data0[1, 2:ncol(data0)])
  datatypes <- c(data0[2, 2:ncol(data0)])
  start <- as.integer(substr(data0[4, 1], 0, 4))
  data <- data.frame(sapply(2:ncol(data0), function(j) as.numeric(unlist(data0[4:nrow(data0), j]))))
  colnames(data) <- colnames(data0)[2:ncol(data0)]

  # nominal to real
  if (makeReal) {
    us <- which(datatypes == "USD")
    for (a in us) {
      data[, a] <- data[, a] / data$CPI_US
    }
  }

  return(list(
    data = data,
    start = start,
    frequency = 12,
    desc = descriptions,
    types = datatypes,
    isReal = makeReal
  ))
}
