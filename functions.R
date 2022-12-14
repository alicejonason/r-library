toeplitz_matrix <- function(x) { # Vector of values used to construct the matrix
  if (length(x) %% 2 == 0) stop("x not of odd length")
  n <- (length(x) + 1) / 2
  toeplitz_mat <- matrix(nrow = n, ncol = n)
  x <- rev(x)
  x <- c(x[-n:-length(x)], rev(x[n:length(x)]))
  for (i in 1:nrow(toeplitz_mat)) {
    for (j in 1:ncol(toeplitz_mat)) {
      toeplitz_mat[i, j] <- x[n - (i - j)]
    }
  }
    return(toeplitz_mat)
}
