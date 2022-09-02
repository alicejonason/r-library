newton_raphson <- function(start_values, # vector of starting values
                           gradient, # function for gradient of f(x_vec)
                           hessian, # function for hessian of f(x_vec)
                           eps = 1e-6) {
  i    <- 0
  cc   <- eps + 100
  x0   <- start_values

  while(cc > eps && i < 1000) {
    x1 <- x0 - solve(hessian(x0)) %*% gradient(x0)
    cc <- sum(abs(x1 - x0)) / (sum(abs(x0 - 0)))
    x0 <- x1
    i  <- i + 1
  }
  return(list(iterations = i, maximum = c(x1)))
}