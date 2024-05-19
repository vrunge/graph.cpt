

p <- 10
n <- 100

data <- matrix(rnorm(n*p),p,n)

ts_to_SPDts(data, 10)


################################################################################


################################################################################

s1 <- 1*sample(dim)
s2 <- 1.5*sample(dim)
s3 <- 2*sample(dim)

s <- matrix(c(s1,s2,s3),dim,3)

n <- 200
node <- c(rep(1, 40), rep(2, 40), rep(3, 20),rep(1, 40), rep(2, 40), rep(3, 20))

res <- NULL

for(i in 1:n)
{
  res <- abind(res,
               Wmatrix(s[,node[i]]),
               along = 3)

}
list_of_matrices <- lapply(seq(dim(res)[3]), function(i) res[, , i])

sigma1 <- diag(s1, dim)
sigma2 <- diag(s2, dim)
sigma3 <- diag(s3, dim)

ll <- list(sigma1,sigma2,sigma3)
ll
distanc <- SPDts_to_dists(list_of_matrices, ll)
dim(distanc)
apply(distanc, 2, which.min)


ei <- eigen(list_of_matrices[[1]])
plot(ei$values, type = 'b')
ei$values



cp <- graph_cpt_manifold(distanc, transition_matrix(3), beta = 0)
plot(cp$best_path, type = 'b', lwd = 2)
par(new = TRUE)
cp <- graph_cpt_manifold(distanc, transition_matrix(3), beta = 0.5)
plot(cp$best_path, type = 'b', col = 3, lwd = 2)
