
set.seed(340)
data <- matrix(rnorm(500), nrow = 50, ncol = 10)
# using prcomp function
resR = prcomp(data, center = TRUE, scale. = TRUE)
# using s.pca in this package
res = s.pca(data,TRUE,TRUE,data)
# res$projections and resR$x must be equal
# res$directions and t(resR$rotation) must be equal

# ----- ANOTHER EXAMPLE: PCA where there is a constant variable:
data <- data.frame( x = rnorm(100), y = rnorm(100), z = rep(0, 100))
res <- s.pca(data)
res_invalid <- try(prcomp(data, center = TRUE,
                          scale. = TRUE)) # Fails, we should remove 'z' first

