# ---------------------------------------------------------------------------- #
# Integral algorithms
# ---------------------------------------------------------------------------- #

# Riemann integral
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
  while (cc > eps) {
    n[t] <- n[t - 1] * 2
    h <- (b - a) / n[t]
    x_i <- a + (0:(n[t] - 1)) * h
    int_apr[t] <- h * sum(f(x_i))
    cc <- abs(int_apr[t] / int_apr[t - 1] - 1)
    t <- t + 1
  }
  return(list(result = tail(int_apr, n = 1),
              nodes = n,
              approximations = int_apr))
}

# Trapezoidal rule
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
  while (cc > eps && t < 2000) {
    n[t] <- n[t - 1] * 2
    h <- (b - a) / n[t]
    x_i <- a + (1:(n[t] - 1)) * h
    int_apr[t] <- h * (f(a) / 2 + sum(f(x_i)) + f(b) / 2)
    cc <- abs(int_apr[t] / int_apr[t - 1] - 1)
    t <- t + 1
  }
  return(list(result = tail(int_apr, n = 1),
              nodes = n,
              approximations = int_apr))
}

# Simpson's rule
simpsons <- function(f, # function
                     a, # lower bound of integral
                     b, # upper bound of integral
                     n) { # number of nodes
  h <- (b - a) / n
  int_apr <- c(0)
  for (i in 1:(n / 2)) {
    int_apr[i] <-
    f(a + h * (2 * i - 2)) + 4 * f(a + h * (2 * i - 1)) + f(a + 2 * h * i)
    }
  return((h / 3) * sum(int_apr))
}

# ---------------------------------------------------------------------------- #
# Root-finding algorithms
# ---------------------------------------------------------------------------- #

# Bisection with absolute convergence criterion
bisection <- function(deriv_fun, # function for computing derivative
                      a,         # lower endpoint
                      b,         # upper endpoint
                      eps = 1e-6) {
  x <- c()
  t <- 1

  if ((deriv_fun(a) * deriv_fun(b)) > 0)
    stop("[a, b] must be s.t. g'(a) * g'(b) < 0")
  x[t] <- (a + b) / 2
  cc  <- eps + 100
  while (cc > eps && t < 100) {
    if ((deriv_fun(a) * deriv_fun(x[t])) <= 0)
      b <- x[t]
    if ((deriv_fun(x[t]) * deriv_fun(b)) < 0)
      a <- x[t]
    x[t + 1] <- (a + b) / 2
    cc <- abs(x[t + 1] - x[t])
    t <- t + 1
  }
  return(list(iterations = t - 1, x_max = tail(x, n = 1)))
}

# Secant method using absolute convergence criterion
secant <- function(deriv_fun, # function for computing derivative
                   x0,        # lower startpoint
                   x1,        # upper startpoint
                   eps = 1e-6) {
  t <- 2
  x <- c()
  x[t - 1] <- x0
  x[t] <- x1
  cc <- eps + 100
  while (cc > eps && t < 100) {
    x[t + 1] <- x[t] - deriv_fun(x[t]) * ((x[t] - x[t - 1]) / 
                                       (deriv_fun(x[t]) - deriv_fun(x[t - 1])))
    cc <- abs(x[t + 1] - x[t])
    t <- t + 1
  }
  return(list(iterations = t - 2, x_max = tail(x, n = 1)))
}

# Newton-Raphson
newton_raphson <- function(start_values, # vector of starting values
                           gradient, # function for gradient of f(x_vec)
                           hessian, # function for hessian of f(x_vec)
                           eps = 1e-6) {
  i    <- 0
  cc   <- eps + 100
  x0   <- start_values

  while (cc > eps && i < 1000) {
    x1 <- x0 - solve(hessian(x0)) %*% gradient(x0)
    cc <- sum(abs(x1 - x0)) / (sum(abs(x0 - 0)))
    x0 <- x1
    i  <- i + 1
  }
  return(list(iterations = i, maximum = c(x1)))
}