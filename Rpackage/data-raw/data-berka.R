## code to prepare `data.berka` dataset

# Specify the location where the raw data is saved:
RawDataFolder <- "D:/Data/Berka/"



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
#' Use this function to load and combine tables from the ‘Berka’ dataset to create a unified data table.
#'
#' @param dirPath A character string representing the path to the downloaded data directory.
#' @param positive A character value determining the positive class. There are four types of loans:
#' 'A' stands for contract finished, no problems, 'B' stands for contract finished,
#' loan not paid, 'C' stands for running contract, OK so far, 'D' stands
#' for running contract, client in debt.
#' @param negative Similar to \code{positive} for negative class.
#' @param rateFun A function to calculate interest rate in loans with the following arguments: \code{amount}, \code{duration}, \code{paymentPerMonth}.
#'
#' @return A data.frame with the following columns:
#' \item{loan_id}{record identifier}
#' \item{status}{original status of the data (A, B, C, or D)}
#' \item{label}{status of paying off the loan transformed to numeric (0,1) by using \code{positive} and \code{negative} arguments. value=1 means default.}
#' \item{amount}{amount of money}
#' \item{payments}{monthly payments}
#' \item{rate}{rates calculated by the \code{rateFun} function}
#' \item{duration_# (#=12,24,36,48,60)}{dummy variables for the duration of the loan}
#' \item{account_frequency_?}{dummy variables for the frequency of issuance of statements. ?="POPLATEK MESICNE" stands for monthly issuance, ?="POPLATEK TYDNE" stands for weekly issuance, ?="POPLATEK PO OBRATU" stands for issuance after transaction.}
#' \item{order_num}{number of payment orders issued for the account of the loan.}
#' \item{order_sum_amount}{sum of amounts of payment orders issued for the account of the loan.}
#' \item{order_related_account_num}{unique number of 'account of the recipient' in payment orders issued for the account of the loan.}
#' \item{order_related_bank_num}{unique number of 'bank of the recipient' in payment orders issued for the account of the loan.}
#' \item{order_has_?}{dummy variables for 'characterization of the payment' in payment orders issued for the account of the loan.}
#' \item{trans_?num}{number of transactions dealt with the account of the loan (in different groups).}
#' \item{trans_?amount_mean}{mean of 'amount of money' in transactions dealt with the account of the loan (in different groups).}
#' \item{trans_?amount_div_balance}{mean of 'amount of money'/'balance after transaction' in transactions dealt with the account of the loan (in different groups).}
#' \item{trans_related_account_num}{unique number of 'account of the partner' in transactions dealt with the account of the loan.}
#' \item{dist_inhabitants_num}{number of inhabitants in the location of the branch of the account of the loan.}
#' \item{dist_muni_#1#2}{number of municipalities with inhabitants #1-#2 in the location of the branch of the account of the loan.}
#' \item{dist_cities_num}{number of cities in the location of the branch of the account of the loan.}
#' \item{dist_ratio_urban_inhabitants}{ratio of urban inhabitants in the location of the branch of the account of the loan.}
#' \item{dist_avg_salary}{average salary in the location of the branch of the account of the loan.}
#' \item{dist_unemployment95}{unemployment rate '95 in the location of the branch of the account of the loan.}
#' \item{dist_unemployment96}{unemployment rate '96 in the location of the branch of the account of the loan.}
#' \item{dist_entrepreneurs_num_per1000}{number of entrepreneurs per 1000 inhabitants in the location of the branch of the account of the loan.}
#' \item{dist_crimes95_num}{number of committed crimes '95 in the location of the branch of the account of the loan.}
#' \item{dist_crimes96_num}{number of committed crimes '96 in the location of the branch of the account of the loan.}
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


# Call the function:
berka <- as.matrix(data.berka.loan(RawDataFolder,
                                          positive = c("B", "D"),
                                          negative = c("A", "C")))

# You might want to use another definition of positive and negative observations:
#data <- Data_BerkaLoan(positive = c("B"), negative = c("A"))

# The dependent variable:
berka_y <- berka[, c("label"), drop = FALSE]
berka_y <- matrix(as.numeric(berka_y), ncol = ncol(berka_y), dimnames = list(NULL, colnames(berka_y)))

# The potential regressors (The first 3 columns of `data` are `id` and `status`, `label`):
berka_x <- berka[, 4:ncol(berka), drop = FALSE]
berka_x <- matrix(as.numeric(berka_x), ncol = ncol(berka_x),
                  dimnames = list(NULL, colnames(berka_x)))
col_means <- colMeans(berka_x, na.rm = TRUE)
large_cols <- col_means > 1000 # These are amounts
berka_x[, large_cols] <- log(berka_x[, large_cols])
colnames(berka_x)[large_cols] <- paste0("ln(", colnames(berka_x)[large_cols],")")

# We use the following weight vector to balance the data
berka_w <- as.matrix((berka_y == 1) * (nrow(berka_y) / sum(berka_y == 1)) + (berka_y == 0))


#descriptions from the documentation:
descriptions <- list(
  list(id="loan_id", description="record identifier"),
  list(id="status", description="original status of the data (A, B, C, or D)"),
  list(id="label", description="default=1"),

  list(id="ln(amount)", description="amount of money (log)"),
  list(id="ln(payments)", description="monthly payments (log)"),
  list(id="rate", description="interest rate"),

  list(id="duration_12", description="dummy variable indicating a 12-month loan duration"),
  list(id="duration_24", description="dummy variable indicating a 24-month loan duration"),
  list(id="duration_36", description="dummy variable indicating a 36-month loan duration"),
  list(id="duration_48", description="dummy variable indicating a 48-month loan duration"),
  list(id="duration_60", description="dummy variable indicating a 60-month loan duration"),

  list(id="account_frequency_POPLATEK_MESICNE", description="dummy variable indicating monthly issuance of statements"),
  list(id="account_frequency_POPLATEK_TYDNE", description="dummy variable indicating weekly issuance of statements"),
  list(id="account_frequency_POPLATEK_PO_OBRATU", description="dummy variable indicating 'after transaction' issuance of statements"),

  list(id="order_num", description="loan account payment orders: number of payment"),
  list(id="ln(order_sum_amount)", description="sum of amounts of payments (log)"),
  list(id="order_related_account_num", description="loan account payment orders: unique number of 'account of the recipient'"),
  list(id="order_related_bank_num", description="loan account payment orders: unique number of 'bank of the recipient'"),

  list(id="order_has_insurrance", description="has insurrance"),
  list(id="order_has_household", description="has household"),
  list(id="order_has_loan", description="has loan"),
  list(id="order_has_leasing", description="has leasing"),



  list(id="trans_?amount_mean", description="mean of 'amount of money' in transactions dealt with the account of the loan (in different groups)."),
  list(id="trans_?amount_div_balance", description="mean of 'amount of money'/'balance after transaction' in transactions dealt with the account of the loan (in different groups)."),
  list(id="trans_related_account_num", description="unique number of 'account of the partner' in transactions dealt with the account of the loan."),
  list(id="dist_inhabitants_num", description="number of inhabitants in the location of the branch of the account of the loan."),

  list(id="dist_muni_0_499", description="number of municipalities in the branch location of the loan account with 1-400 inhabitants"),
  list(id="dist_muni_#1#2", description="number of municipalities in the branch location of the loan account with #1-#2 inhabitants"),


  list(id="dist_cities_num", description="number of cities in the location of the branch of the account of the loan."),
  list(id="dist_ratio_urban_inhabitants", description="ratio of urban inhabitants in the location of the branch of the account of the loan."),
  list(id="dist_avg_salary", description="average salary in the location of the branch of the account of the loan."),
  list(id="dist_unemployment95", description="unemployment rate '95 in the location of the branch of the account of the loan."),
  list(id="dist_unemployment96", description="unemployment rate '96 in the location of the branch of the account of the loan."),
  list(id="dist_entrepreneurs_num_per1000", description="number of entrepreneurs per 1000 inhabitants in the location of the branch of the account of the loan."),
  list(id="dist_crimes95_num", description="number of committed crimes '95 in the location of the branch of the account of the loan."),
  list(id="dist_crimes96_num", description="number of committed crimes '96 in the location of the branch of the account of the loan.")
)


data.berka = list(y = berka_y, x = berka_x, w = berka_w, descriptions = descriptions)


usethis::use_data(data.berka, overwrite = TRUE)
