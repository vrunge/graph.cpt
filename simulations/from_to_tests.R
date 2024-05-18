

p <- 10
n <- 100

data <- matrix(rnorm(n*p),p,n)

ts_to_SPDts(data, 10)
