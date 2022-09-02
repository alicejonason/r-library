trapezoidal <- function(f, # function
                        a, # lower bound of integral
                        b, # upper bound of integral
                        n, # number of nodes
                        eps = 1e-6) {
  n <- c(n)
  h <- (b - a) / n
  x_i <- a + (1:(n - 1)) * h
  int_apr <- c(h * (f(a) / 2 + sum(f(x_i)) + f(b) / 2))
  t <- 2
  cc <- eps + 100
  while(cc > eps && t < 2000) {
    n[t] <- n[t - 1] * 2
    h <- (b - a) / n[t]
    x_i <- a + (1:(n[t] - 1)) * h
    int_apr[t] <- h * (f(a) / 2 + sum(f(x_i)) + f(b) / 2)
    cc <- abs(int_apr[t] / int_apr[t - 1] - 1)
    t <- t + 1
  }
  return(list(result = tail(int_apr, n = 1), nodes = n, approximations = int_apr))
}
