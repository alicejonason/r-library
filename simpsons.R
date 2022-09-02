simpsons <- function(f, # function
                     a, # lower bound of integral
                     b, # upper bound of integral
                     n) { # number of nodes#
  h <- (b - a) / n
  int_apr <- c(0)
  for (i in 1:(n / 2)) {
    int_apr[i] <-
    f(a + h * (2 * i - 2)) + 4 * f(a + h * (2 * i - 1)) + f(a + 2 * h * i)
    }
  return((h / 3) * sum(int_apr))
}
