riemann <- function(f, # function
                    a, # lower bound of integral
                    b, # upper bound of integral
                    n, # number of nodes
                    eps = 1e-6) {
  n <- c(n)
  h <- (b - a) / n
  x_i <- a + (0:(n - 1)) * h
  int_apr <- c(h * sum(f(x_i)))
  t <- 2
  cc <- eps + 100
  while(cc > eps) {
    n[t] <- n[t - 1] * 2
    h <- (b - a) / n[t]
    x_i <- a + (0:(n[t] - 1)) * h
    int_apr[t] <- h * sum(f(x_i))
    cc <- abs(int_apr[t] / int_apr[t - 1] - 1)
    t <- t + 1
  }
  return(list(nodes = n, approximations = int_apr))
}
