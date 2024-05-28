
vec <- c(1,2,3)
D <- diag(vec)
n <- 3

random_matrix <- matrix(rnorm(n^2), n, n)
Q <- qr.Q(qr(random_matrix)) # Q%*%t(Q) = Id
t(Q)%*%Q
solve(Q)%*%Q
W <- Q %*% D %*% t(Q)

eigen(W)
eigen(solve(W))


### ### ### ### ### SAME RESULTS ### ### ### ### ###

dist_SPD(W, diag(vec))
dist_SPD(solve(W), solve(diag(vec)))
dist_SPD(Q %*% W %*% t(Q), Q %*% diag(vec) %*% t(Q))


