## Code for preparing `data.pcp` dataset
## This uses 'IMF's Primary Commodity Prices' dataset
## You need to download the Excel file from the website
## The frequency of data is monthly

filePath <- "D:/Data/PCP/external-datasep.xls"

data0 <- readxl::read_excel(filePath) # update this

descriptions <- c(data0[1, 2:ncol(data0)])
datatypes <- c(data0[2, 2:ncol(data0)])
start <- as.integer(substr(data0[4, 1], 0, 4))

data.pcp <- data.frame(sapply(2:ncol(data0), function(j) as.numeric(unlist(data0[4:nrow(data0), j]))))
colnames(data.pcp) <- colnames(data0)[2:ncol(data0)]
rownames(data.pcp) <- tdata::get.seq0(tdata::f.monthly(start,1),nrow(data.pcp))

# Just save the Indices:
n <- 18
data.pcp <- list(data = data.pcp[,1:n],
                 descriptions = descriptions[1:n],
                 datatypes = unlist(datatypes)[1:n],
                 start = start)

usethis::use_data(data.pcp, overwrite = TRUE)
