dirichlet_sim <- function(n, alpha) {
  theta <- as.data.frame(matrix(NA, n, length(alpha)))
  for (j in 1:length(alpha)) {
    theta[, j] <- rgamma(n, alpha[j], 1)
  }
  for (i in 1:n) {
    theta[i, ] <- theta[i, ] / sum(theta[i, ])
  }
  return(theta)
}