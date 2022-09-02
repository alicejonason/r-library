dirichlet_sim <- function(n_iter, alpha) {
  n_cat <- length(alpha)
  theta <- matrix(NA, n_iter, n_cat)
  for (j in 1:n_cat) {
    theta[, j] <- rgamma(n_iter, alpha[j], 1)
  }
  for (i in 1:n_iter) {
    theta[i, ] <- theta[i, ] / sum(theta[i, ])
  }
  return(theta)
}