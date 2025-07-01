## period 2: sample baseline hazard rate parameter
lambda1_Gibbs <- function(n.censor.t1, time.1j){
  a <- n.censor.t1 + 1
  b <- time.1j + 1
  lambda1 <- rgamma(1, shape = a, rate = b)
  return(lambda1)
}

## period 2
# sample baseline hazard rate parameter
lambda_Gibbs <- function(n.censor.all, time.j, u){
  a <- sum(n.censor.all) + 1
  b <- time.j$time
  b <- sum(b * u) + 1
  lambda <- rgamma(n = 1, shape = a, rate = b)
  return(lambda)
}

# sample center specific random effect uj here
u_Gibbs <- function(n.censor, time.j, lambda, theta){
  a <- n.censor$ncensor + theta
  Tj <- time.j$time
  b <- lambda * Tj + theta
  params <- matrix(c(a, b), ncol = 2, byrow = FALSE)
  u <- apply(params, 1, FUN = function(x){rgamma(n = 1, shape = x[1], rate = x[2])})
  return(u)
}

# sample center specific random probability vector p here
p_Gibbs <- function(n.j, gamma){
  params <- as.matrix(n.j)
  params <- cbind(params, gamma)
  p <- apply(params, 1, FUN = function(x){rdirichlet(n = 1, alpha = x[1:3] + x[4:6])})
  return(t(p))
}

