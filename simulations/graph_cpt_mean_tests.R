

vec <- c(rnorm(10), rnorm(10, mean = 5), rnorm(10, mean = 8))
vec

res <- graph_cpt_mean(y = vec,
                      A = transition_matrix(3),
                      states = c(0,5,10),
                      beta = 0)
res
length(res[[3]])

