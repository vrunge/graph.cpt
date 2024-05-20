
vec <- c(1,2,3)
D <- diag(vec)
n <- 3
random_matrix <- matrix(rnorm(n^2), n, n)
Q <- qr.Q(qr(random_matrix)) # Q%*%t(Q) = Id
t(Q)%*%Q
solve(Q)%*%Q
W <- Q %*% D %*% t(Q)
eigen(W)
eigen(diag(vec))


cat("res" ,dist_SPD(D,diag(vec)))
cat("res : " ,dist_SPD(W,diag(vec)), " ")
cat("res : " ,dist_SPD(W,diag(c(vec[3],vec[2],vec[1]))), " ")



