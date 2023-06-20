sample1 <- sim.varma(2L, 3L, 2L)

ar1 <- matrix(c(0.7,0.2,-0.4,0.3),2,2)
ar2 <- matrix(c(-0.4,0.1,0.2,-0.3),2,2)
ma1 <- matrix(c(0.5,-0.1,0.3,0.4),2,2)
Sigma <- matrix(c(1,0.3,0.3,1),2,2)
B <- matrix(c(0.5,-0.3),2)

sample2 <- sim.varma(Sigma, list(ar1, ar2), list(ma1), exoCoef = B ,
                    nObs =100, nBurn =10 , intercept = c(1,-1))

# Plot the y series
matplot(sample2$y,type = "l")
