
vec <- c(1,2,3)
D <- diag(vec)
n <- 3
random_matrix <- matrix(rnorm(n^2), n, n)

Q <- qr.Q(qr(random_matrix)) # Q%*%t(Q) = Id

t(Q)%*%Q
solve(Q)%*%Q
W <- Q %*% D %*% t(Q)




library(rWishart)
data <- rWishart(10, 10, diag(1, 5))
for(i in 1:10) print(dist_SPD(data[,,1], data[,,i]))

data <- rWishart(10, df = 10^5, diag(1, 5))
for(i in 1:10) print(dist_SPD(data[,,1], data[,,i]))

data <- rWishart(10, 10^10, diag(1, 5))
for(i in 1:10) print(dist_SPD(data[,,1], data[,,i]))


