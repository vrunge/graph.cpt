
library(expm)
n <- 5
p <- 10
M1 <- matrix(rnorm(p*n), n, p)
S1 <- M1%*%t(M1)
eigen(S1)$values

S <- solve(sqrtm(S1))%*%S1%*%solve(sqrtm(S1))
sqrt(sum(log(eigen(S)$values)^2))
dist_SPD(S1,S1)


###

M2 <- matrix(rnorm(p*n), n, p)
S2 <- M2%*%t(M2)
S <- solve(sqrtm(S1))%*%S2%*%solve(sqrtm(S1))
sqrt(sum(log(eigen(S)$values)^2))
dist_SPD(S1,S2)


library(expm)
n <- 5
p <- 10
M1 <- matrix(rnorm(p*n), n, p)
M2 <- matrix(rnorm(p*n), n, p)
dist_SPD(S1,S2)











