# Function to calculate the derivative of a function x
derivative <- function(f, # function to take derivative of
                       var_name) { # variable (char) to take derivative w.r.t.
  D(parse(text = as.character(body(f)[2])), var_name)
}