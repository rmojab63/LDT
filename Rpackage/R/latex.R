
get_coef_stars <- function(pvalue, formatLatex) {
  if (is.nan(pvalue)) { # e.g., restricted to zero
    return(if (formatLatex) "\\textsuperscript{(r)}" else "<sup>(r)</sup>")
  }
  paste0(
    (if (formatLatex) "\\textsuperscript{" else "<sup>"),
    if (pvalue <= 0.01) {
      "***"
    } else if (pvalue <= 0.05) {
      "**"
    } else if (pvalue <= 0.1) {
      "*"
    } else {
      ""
    }, (if (formatLatex) "}" else "</sup>")
  )
}

get_coef_func <- function(tableFun, formatNumFun, formatLatex){

  if (is.function(tableFun)){
    # no action
  }
  else if (tableFun == "sign") {
    tableFun <- function(j, coef, std, pvalue) {
      if (coef > 0) {
        "+"
      } else if (coef < 0) {
        "-"
      } else {
        "0"
      }
    }
  }
  else if (tableFun == "sign_star") {
    tableFun <- function(j, coef, std, pvalue) {
      paste0(if (coef > 0) {
        "+"
      } else if (coef < 0) {
        "-"
      } else {
        "0"
      }, get_coef_stars(pvalue, formatLatex))
    }
  }
  else if (tableFun == "coef") {
    tableFun <- function(j, coef, std, pvalue) {
      formatNumFun(j, coef)
    }
  }
  else if (tableFun == "coef_star") {
    tableFun <- function(j, coef, std, pvalue) {
      paste0(formatNumFun(j, coef), get_coef_stars(pvalue, formatLatex))
    }
  }
  else if (tableFun == "coef_star_std") {
    tableFun <- function(j, coef, std, pvalue) {
      paste0(formatNumFun(j, coef), get_coef_stars(pvalue, formatLatex), " (", formatNumFun(j, std), ")")
    }
  }
  else
    stop("tableFun must be a function or a valid character string.")

  tableFun
}

#' Create Table of Coefficients
#'
#' This function summarizes a list of estimated models (output of \code{estim.?} functions) and creates
#' a table of coefficients.
#'
#' @param estimList A named list where each element is output from a \code{estim.?} function, all belonging to a common analysis.
#' @param depList List of dependent variable name to be included in the columns of the table. If \code{NULL}, everything is added.
#' @param tableFun Function with arguments \code{coef}, \code{std}, \code{pvalue}, \code{minInColm}, \code{maxInCol} for formatting estimated coefficients.
#' Can also use one of the following character strings for predefined formatting: \code{sign}, \code{sign_star}, \code{coef}, \code{coef_star}, \code{coef_star_std}.
#' @param formatNumFun Function to format numbers if \code{tableFun} uses predefined formatting. It takes two arguments: \code{colIndex} (determines the column index) and \code{x} (determines the value).
#' @param regInfo List of keys (such as \code{num_eq}, \code{num_dep}, ...) to determine information at bottom of table. Use "" (empty) for empty rows.
#' A list of available options is given in the details section.
#' @param textFun Function to change any text in columns or rows of the table to a more informative text.
#' It has two arguments: \code{text} and \code{type}.
#' @param textFun_sub List for replacing special characters. If \code{NULL}, 'list(c("%", "\\\\%"), c("_", "\\\\_"))' is used.
#' @param textFun_max Maximum length for texts in the table.
#' @param expList Determines the name of the explanatory variables to insert in table.
#' If \code{NA}, it inserts all coefficients.
#' If it is an integer, it insert that number of explanatory variables in the table.
#' It can be a list of available explanatory variables.
#' Use \code{...} for empty rows.
#' @param formatLatex If \code{TRUE}, default options are for 'latex', otherwise, 'html'.
#' @param numFormat default formatting for the numbers.
#'
#' @details
#' The first part of the table is the header, followed by the coefficients. At the bottom, you can insert
#' the following items by specifying \code{regInfo}:
#' \itemize{
#' \item An empty character string (i.e., "") for inserting empty line.
#' \item \code{"sigma2"} for the covariance of regression, if it is available.
#' \item Name of an element in \code{estimList[[...]]$counts}.
#' \item An available metric name in the row names of \code{estimList[[...]]$metrics}.
#' }
#'
#' Furthermore, second argument in \code{textFun} can be:
#' \itemize{
#' \item \code{hname}: shows that the text is a header name from the \code{estimList} elements.
#' \item \code{dname}: shows that the text is a dependent variable name from the columns of \code{coefs} matrix.
#' \item \code{rname}: shows that the text is a key given in \code{regInfo}.
#' \item \code{ename}: shows that the text is an explanatory variable name from the rows of \code{coefs} matrix.
#' \item \code{NULL}: shows that the text is a specific code or something else.
#' }
#'
#' @return A data frame with (formatted) requested information.
#' @export
#' @examples
#' # See 'search.?' or 'estim.?' functions for some examples.
#'
coefs.table <- function(estimList, depList = NULL, tableFun = "coef_star", formatNumFun  = NULL,
                        regInfo = NULL, textFun = NULL,
                        textFun_sub = NULL, textFun_max = 20,
                        expList = NA, formatLatex = TRUE,
                        numFormat = "%.2f") {

  if (is.null(estimList$counts) == FALSE)
    estimList <- list(m = estimList)


  if (is.null(regInfo))
    regInfo <- c("", "obs")

  # there are some functions that format the inputs:
  # Function to change or truncate the text:
  if (is.null(textFun))
    textFun <- function(text, type) text
  textFunCopy = textFun
  textFun <- function(text, type){
    text <- textFunCopy(text, type)
    if (nchar(text) > textFun_max)
      text <- paste0(substr(text, 1, textFun_max - 3), "...")
    text
  }

  # Function to format the numbers

  if (is.null(formatNumFun))
    formatNumFun = function(colIndex, x) {
      if (is.integer(x))
        x
      else
        sprintf0(numFormat, x)
    }


  # Function to replace special characters
  if(is.null(textFun_sub))
    textFun_sub = list(c("%", "\\\\%"), c("_", "\\\\_"))
  textFun_sub <- as.list(textFun_sub)
  textSubFun = function(text){
    for (sg in textFun_sub)
      text <- gsub(sg[[1]], sg[[2]], text, fixed = TRUE)
    text
  }

  # Function to format the coefficients:
  tableFun <- get_coef_func(tableFun, formatNumFun, formatLatex)

  # Create column names
  col_names <- names(estimList)
  if (is.null(col_names))
    col_names <- paste0("m", rep(1:length(estimList)))
  col_names <- lapply(col_names, function(n) textFun(n, "hname"))

  if (is.null(depList)) # get all names
    depList <- unique(unlist(lapply(estimList, function(e) colnames(e$estimations$coefs))))
  else
    depList <- as.character(depList)
  depList <- lapply(depList, function(n) textFun(n, "dname"))

  col_names <- unlist(lapply(seq_along(col_names), function(i) {
    rep(col_names[i], sum(colnames(estimList[[i]]$estimations$coefs) %in% depList))
  }))

  # Create row names
  if (is.character(expList))
    row_names <- expList
  else{
    row_names <- unique(unlist(lapply(estimList, function(e) row.names(e$estimations$coefs)))) # get all names
    if (length(expList) == 1 && is.numeric(expList))
      row_names <- unlist(c(row_names[1:expList], "..."))
  }
  row_names0 <- row_names # without formatting
  last_x_ind <- length(row_names)
  #    include other rows:
  dep_cel_code <- "dep."
  row_names <- unlist(c(textFun(dep_cel_code, NULL),
                        lapply(row_names, function(n) textFun(n, "ename")),
                        lapply(regInfo, function(n) textFun(n, "rname"))))
  row_names0 <- unlist(c(dep_cel_code, row_names0, regInfo))


  # Create the table
  r_table <- matrix("", length(row_names), length(col_names))
  rownames(r_table) <- row_names
  colnames(r_table) <- col_names

  # start the loop
  j <- 0
  for (e in estimList) {
    for (d in depList){
      if (!d %in% colnames(e$estimations$coefs)) next
      j <- j + 1
      r_table[1, j] <- d

      ns <- row.names(e$estimations$coefs)
      i <- 0

      # Insert coefficients:
      for (v in ns) {
        i <- i + 1

        ind <- which(row_names0 == v)
        if (length(ind) == 0) # not found in the table
          next
        if (ind[[1]] > last_x_ind+1) # it's a regInfo name
          next
        if (length(ind) > 1)
          warning(paste0("Multiple name exists (name= ", v ,"estim. index = ", j, ")"))

        if (v == "...") # allow it
          next

        coef <- e$estimations$coefs[i, d]
        std <- e$estimations$stds[i, d]
        pv <- e$estimations$pValues[i, d]
        r <- tableFun(j, coef, std, pv)
        r_table[ind[[1]], j] <- r
      }

      # Insert regression information:
      i <- last_x_ind + 1 # because of the first row
      for (rn in regInfo) {
        i <- i + 1
        v <- NULL
        r <- rn[[1]]

        if (r == "")
          v <- "" # for an empty cell/line
        else if (r == "sigma2")
          v <- formatNumFun(j, e$estimations$sigma[d, d])
        else{
          # check and see if it is a metric
          ind <- which(r == rownames(e$metrics))
          if (length(ind) != 0) {
            v <- formatNumFun(j, e$metrics[ind, d])

            if (r == "f") {
              ind0 <- which("fProb" == rownames(e$metrics))
              fp <- e$metrics[ind0, d]
              if (is.na(fp)) {
                warning("p-value of 'F' statistics is NA.")
              }
              v <- paste0(v, get_coef_stars(fp))
            }

          }
          else{

            # check and see if it is a count metric
            ind <- which(r == names(e$counts))
            if (length(ind) != 0) {
              v <- formatNumFun(j, e$counts[[ind]])
            }

          }
        }

        if (is.null(v) || length(v) == 0) {
          warning(paste0("Invalid 'regInfo' member is requested. 'NA' is used. code=", r))
          v <- NA
        }
        if (ncol(r_table) < j) {
          d <- 0
        }

        r_table[i, j] <- v

      }
    }
  }

  if (length(estimList) == 1){ # header is the same for all
    colnames(r_table) <- r_table[1,]
    r_table <- r_table[2:nrow(r_table),,drop=FALSE]
  }
  else if (length(depList) == 1){
    r_table <- r_table[2:nrow(r_table),,drop=FALSE]
  }

  return(r_table)
}


sprintf0 <- function(numFormat, x) {
  if (length(x) == 0) {
    return(ifelse(x==0,"0",sprintf(numFormat, x)))
  } else {
    return(sapply(x, function(d)ifelse(d==0,0,sprintf(numFormat, d))))
  }
}


#' Convert a matrix to LaTeX code
#'
#' This function takes a matrix and returns the LaTeX code for the matrix,
#' using the specified matrix environment and number format.
#'
#' @param mat A matrix to convert to LaTeX code.
#' @param env A character string specifying the matrix environment to use.
#'   Possible values are "bmatrix", "pmatrix", "Bmatrix", "vmatrix", and "Vmatrix".
#' @param numFormat A character string specifying the format to use when displaying the numbers.
#'   This value can be any valid format string accepted by the sprintf function.
#'
#' @return A character string containing the LaTeX code for the matrix.
#'
#' @examples
#' #mat <- matrix(1:4, ncol = 2)
#' #latex.matrix(mat)
#' #latex.matrix(mat, env = "pmatrix")
#' #latex.matrix(mat, numFormat = "%.0f")
latex.matrix <- function(mat, env = "bmatrix", numFormat = "%.2f") {

  mat_str <- apply(mat, 1, function(row) {
    paste0(sprintf0(numFormat, row), collapse = " & ")
  })
  mat_str <- paste(mat_str, collapse = " \\\\ ")
  mat_str <- paste0("\\begin{", env, "} ", mat_str, " \\end{", env, "}")
  return(mat_str)
}

#' Generate LaTeX code for a variable vector
#'
#'
#' @param vec_size An integer specifying the size of the vector.
#' @param label A character string specifying the label to use when filling the cell values.
#' @param intercept A logical value indicating whether to add `1` to the vector.
#' @param max_size An integer specifying the maximum size of the vector before using `\\vdots` to remove middle elements.
#' @param env A character string specifying the matrix environment to use.
#'   Possible values are "bmatrix", "pmatrix", "Bmatrix", "vmatrix", and "Vmatrix".
#'
#' @return A character string containing the LaTeX code for the vector.
#'
#' @examples
#' #latex.variable.vector(3, "X")
#' #latex.variable.vector(5, "X")
#' #latex.variable.vector(5, "X", intercept = FALSE)
latex.variable.vector <- function(vec_size, label, intercept = FALSE, max_size = 3, env = "pmatrix") {


  if (vec_size == 0){
    if (intercept)
      vec_names <- paste0("\\begin{",env,"}   1 \\end{",env,"}")
    else
      vec_names <- paste0("\\begin{",env,"} \\end{",env,"}")
  }
  else if (vec_size <= max_size) {
    vec_names <- paste0(label,"_", seq_len(vec_size - ifelse(intercept,1,0)), collapse = " \\\\ ")
    if (intercept) {
      vec_names <- paste0("1 \\\\ ", vec_names)
    }
    vec_names <- paste0("\\begin{",env,"} ", vec_names, " \\end{",env,"}")
  } else {
    vec_names <- paste0("\\begin{",env,"}")
    if (intercept) {
      vec_names <- paste0(vec_names,"1 \\\\ ",label,"_1 \\\\ \\vdots \\\\ ",label,"_", vec_size - 1)
    } else {
      vec_names <- paste0(vec_names,label,"_1 \\\\ \\vdots \\\\ ",label,"_", vec_size)
    }
    vec_names <- paste0(vec_names,"\\end{",env,"}")
  }

  return(vec_names)
}



sur.to.latex.eqs <- function(sigma, coef, intercept, numFormat = "%.2f") {
  num_y <- ncol(coef)
  num_x <- nrow(coef) - ifelse(intercept, 1, 0)
  eqs <- character(num_y)
  for (i in seq_len(num_y)) {
    b <- coef[, i]
    if (intercept) {
      eqs[i] <- paste0("Y_", i, " = ", sprintf0(numFormat, b[1]))
      if (num_x > 0) {
        eqs[i] <- paste0(eqs[i], " + ", paste0(sprintf0(numFormat, b[-1]), " X_", seq_len(num_x), collapse = " + "))
      }
    } else {
      eqs[i] <- paste0("Y_", i, " = ", paste0(sprintf0(numFormat, b), " X_", seq_len(num_x), collapse = " + "))
    }
    eqs[i] <- paste0(eqs[i], " + E_", i, ", \\sigma_",i,"^2 = ", sprintf0(numFormat, sigma[[i,i]]))
  }
  eqs_latex <- paste(eqs, collapse = " \\\\ ")

  return(eqs_latex)
}

sur.to.latex.mat <- function(sigma, coef, intercept = TRUE, numFormat = "%.2f",
                             num_x_break = 3, y_label = "Y", x_label= "X", e_label = "E") {
  num_y <- ncol(coef)
  num_x <- nrow(coef)

  y_vec <- latex.variable.vector(num_y, y_label)

  coef_t <- t(coef)
  x_mat <- latex.matrix(mat = coef_t, numFormat = numFormat)

  x_vec <- latex.variable.vector(num_x, x_label, intercept, num_x_break)

  e_vec <- latex.variable.vector(num_y, e_label)

  s_mat <- latex.matrix(mat = sigma, numFormat = numFormat)

  eq_latex <- paste0(y_vec, " = ", x_mat," ", x_vec," + ", e_vec, ", \\Sigma = ", s_mat)



  return(eq_latex)
}




varma.to.latex.mat <- function(sigma, arList, int, exoCoef, maList, d, D, s, numFormat = "%.2f") {

  #TODO: add breaks based on the number of equations and lags
  #TODO: add three vertical dots for large vectors similar to SUR

  numEq <- nrow(sigma)
  numAR <- ifelse(is.null(arList), 0 , length(arList))
  numMA <- ifelse(is.null(maList), 0 , length(maList))
  numExo <- ifelse(is.null(exoCoef), 0 , ncol(exoCoef))

  # Initialize the LaTeX string
  latex_str <- ""

  delta <- ""
  if (d == 1)
    delta <- paste0(delta, "\\Delta")
  else if (d > 1)
    delta <- paste0(delta, "\\Delta^", d)
  if (D == 1)
    delta <- paste0(delta, "\\Delta_", s)
  else if (D > 1)
    delta <- paste0(delta, "\\Delta_", s, "^", D)

  # Add the dependent variable vector
  latex_str <- paste0(latex_str, " \\begin{bmatrix}", paste(paste0(delta, " Y_{", seq_len(numEq),"t}"), collapse = "\\\\"), "\\end{bmatrix}")

  # Add the equal sign
  latex_str <- paste0(latex_str, " = ")

  # Add the intercept vector
  if (!is.null(int) && any(as.numeric(int) != 0)) {
    latex_str <- paste0(latex_str, "\\begin{bmatrix}", paste(sprintf(numFormat, int), collapse = "\\\\"), "\\end{bmatrix}")

    # Add a plus sign if there are more terms
    if (numAR > 0 || numExo > 0 || numMA > 0) {
      latex_str <- paste0(latex_str, " + ")
    }
  }

  # Add the AR matrices
  for (lag in seq_len(numAR)) {
    if (all(as.numeric(arList[[lag]])== 0))
      next
    ar_mat <- matrix(sprintf(numFormat, arList[[lag]]), nrow = nrow(arList[[lag]]), ncol = ncol(arList[[lag]]))
    latex_str <- paste0(
      latex_str,
      "\\begin{bmatrix}",
      paste(apply(ar_mat, 1, paste, collapse = " & "), collapse = "\\\\"),
      "\\end{bmatrix}"
    )

    # Add the lagged dependent variable vector
    latex_str <- paste0(latex_str,
                        " \\begin{bmatrix}",
                        paste(paste0(delta, " Y_{", seq_len(numEq), "t-", lag, "}"), collapse = "\\\\"),
                        "\\end{bmatrix}"
    )

    # Add a plus sign if there are more terms
    if (lag < numAR || numExo > 0 || numMA > 0) {
      latex_str <- paste0(latex_str, " + ") #break
    }
  }

  if (numExo > 1){
    latex_str <- paste0(latex_str, "\\\\") #break
  }

  # Add the exogenous matrix
  if (numExo > 0 && any(as.numeric(exoCoef) != 0)) {
    exo_mat <- matrix(sprintf(numFormat, exoCoef), nrow = nrow(exoCoef), ncol = ncol(exoCoef))
    latex_str <- paste0(
      latex_str,
      "\\begin{bmatrix}",
      paste(apply(exo_mat, 1, paste, collapse = " & "), collapse = "\\\\"),
      "\\end{bmatrix}"
    )

    # Add the exogenous variable vector
    if (numExo > 0){
      latex_str <- paste0(
        latex_str,
        "\\begin{bmatrix}",
        paste(paste0("X_", seq_len(numExo)), collapse = "\\\\"),
        "\\end{bmatrix}"
      )
    }

    # Add a plus sign if there are more terms
    if (numMA > 0) {
      latex_str <- paste0(latex_str, "+ ")
    }
  }

  if (numMA > 0){
    latex_str <- paste0(latex_str, "\\\\") #break
  }

  # Add the MA matrices
  for (lag in seq_len(numMA)) {
    if (all(as.numeric(maList[[lag]])== 0))
      next
    ma_mat <- matrix(sprintf(numFormat, maList[[lag]]), nrow = nrow(maList[[lag]]), ncol = ncol(maList[[lag]]))
    latex_str <- paste0(
      latex_str,
      "\\begin{bmatrix}",
      paste(apply(ma_mat, 1, paste, collapse = " & "), collapse = "\\\\"),
      "\\end{bmatrix}"
    )

    # Add the lagged error vector
    latex_str <- paste0(
      latex_str,
      "\\begin{bmatrix}",
      paste(paste0("E_{", seq_len(numEq), "t-", lag, "}"), collapse = "\\\\"),
      "\\end{bmatrix}"
    )

    # Add a plus sign if there are more terms
    if (lag < numMA) {
      latex_str <- paste0(latex_str, " + ")
    }
  }

  # Add the error vector
  latex_str <- paste0(latex_str, " + \\begin{bmatrix}", paste(paste0("E_{", seq_len(numEq),"t}"), collapse = "\\\\"), "\\end{bmatrix}")

  s_mat <- latex.matrix(mat = sigma, numFormat = numFormat)
  latex_str <- paste0(latex_str, ",\\\\ \\Sigma = ", s_mat)

  return(latex_str)
}

varma.to.latex.eqs <- function(sigma, arList, int, exoCoef, maList, d, D, s, numFormat = "%.2f") {

  #TODO: handle zero coefficients

  numEq <- nrow(sigma)
  numAR <- ifelse(is.null(arList), 0 , length(arList))
  numMA <- ifelse(is.null(maList), 0 , length(maList))
  numExo <- ifelse(is.null(exoCoef), 0 , ncol(exoCoef))

  # Initialize the LaTeX string
  latex_str <- ""

  delta <- ""
  if (d == 1)
    delta <- paste0(delta, "\\Delta")
  else if (d > 1)
    delta <- paste0(delta, "\\Delta^", d)
  if (D == 1)
    delta <- paste0(delta, "\\Delta_", s)
  else if (D > 1)
    delta <- paste0(delta, "\\Delta_", s, "^", D)

  # Add the equations
  for (eq in seq_len(numEq)) {
    # Add the dependent variable
    latex_str <- paste0(latex_str, delta," Y_{", eq, "t}")

    # Add the intercept
    if (!is.null(int) && int[eq] != 0) {
      latex_str <- paste0(latex_str, " = ", sprintf0(numFormat, int[eq]))
    }

    # Add the AR terms
    for (lag in seq_len(numAR)) {
      coefs <- arList[[lag]][eq, ]
      if (all(as.numeric(coefs) == 0))
        next
      signs <- ifelse(coefs >= 0, " + ", " - ")
      coefs <- sprintf0(numFormat, abs(coefs))

      terms <- paste0(signs, coefs, delta, " Y_{", seq_len(numEq), "t-", lag, "}")
      latex_str <- paste0(latex_str, paste(terms, collapse = ""))
    }

    # Add the exogenous terms
    if (numExo > 0){
      for (exo in seq_len(numExo)) {
        coef <- exoCoef[eq, exo]
        if (all(as.numeric(coefs) == 0))
          next
        sign <- ifelse(coef >= 0, " + ", " - ")
        coef <- sprintf0(numFormat, abs(coef))
        term <- paste0(sign, coef, " X_", exo)
        latex_str <- paste0(latex_str, term)
      }
    }

    # Add break
    latex_str <- paste0(latex_str, "\\\\")

    # Add the MA terms
    for (lag in seq_len(numMA)) {
      coefs <- maList[[lag]][eq, ]
      if (all(as.numeric(coefs) == 0))
        next
      signs <- ifelse(coefs >= 0, " + ", " - ")
      coefs <- sprintf0(numFormat, abs(coefs))
      terms <- paste0(signs, coefs, " E_{", seq_len(numEq), "t-", lag, "}")
      latex_str <- paste0(latex_str, paste(terms, collapse = ""))
    }

    # Add the error term
    latex_str <- paste0(latex_str, " + E_{", eq, "t},\\quad \\sigma_", eq,"^2 = ",
                        sprintf0(numFormat, sigma[[eq,eq]]))

    # Add a line break if this is not the last equation
    if (eq < numEq)
      latex_str <- paste0(latex_str, "\\\\")
  }

  return(latex_str)
}



bin.to.latex.eq <- function(coef, probit, numFormat = "%.2f") {

  xNames <- paste0("X_",c(1:length(coef)))

  terms <- character(length(coef))
  terms[1] <- sprintf0(numFormat, coef[1])
  for (i in seq_along(coef)[-1]) {
    # avoid +- signs
    if (coef[i] >= 0) {
      sign <- " + "
    } else {
      sign <- " - "
    }
    terms[i] <- paste0(sign, sprintf(numFormat, abs(coef[i])), " ", xNames[i])
  }

  formula_str <- paste(terms, collapse = "")

  cond_str = ifelse(length(xNames) == 4, paste(xNames[-1], collapse = ", "), paste(c(xNames[2], "...",xNames[length(xNames)]), collapse = ", "))

  if (probit) {
    formula_str <- paste0("P(Y = 1 | ", cond_str, ") = \\Phi(", formula_str, ")")
  } else {
    formula_str <- paste0("P(Y = 1 | ", cond_str, ") = \\frac{1}{1 + e^{-(", formula_str, ")}}")
  }

  return(formula_str)
}
