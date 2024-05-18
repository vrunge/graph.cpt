

vec <- c(rnorm(10, sd = 1), rnorm(10, mean = 1, sd = 0), rnorm(10, mean = 1.9, sd = 1))
vec
states <- c(0,1,2)

res <- graph_cpt_mean(y = vec,
                      A = transition_matrix(3, type = "cyclic"),
                      states = states,
                      beta = 0)
res

plot(vec, type = 'b', ylim = c(min(vec), max(vec)))
par(new = TRUE)
plot(states[res[[3]]], ylim = c(min(vec), max(vec)), col = 2, lwd = 3, type = 'b')


res$path



