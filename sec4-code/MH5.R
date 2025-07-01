# update theta and r together with one AM 
#library(mvtnorm)
#========================================================
# theta
# log-likelihood for all centers
sum_center_loglik_theta <- function(u, theta){
  theta.list <- rep(theta, length(u))
  loglik <- -lgamma(theta.list) + theta.list * log(theta.list) + (theta.list - 1) * log(u) - theta.list * u
  return(sum(loglik))
}
# r
# log-likelihood for each center
center_loglik <- function(pj, r){
  # pj = (p.pos.j, p.neg.j, p0.j)
  # r = (r1, r2, r3)
  p.pos <- pj[1]
  p.neg <- pj[2]
  p0 <- pj[3]
  r1 <- r[1]
  r2 <- r[2]
  r3 <- r[3]
  # check the case when p == 0 ?
  if(sum(pj == 0) > 0){
    print('warning, p == 0')
    stop
  }else{
    loglik <- lgamma(r1 + r2 + r3) - 
      lgamma(r1) - lgamma(r2) - lgamma(r3) + 
      (r1 - 1) * log(p.pos) + (r2 - 1) * log(p.neg) + (r3 - 1) * log(p0)
  }
  return(loglik)
}

sum_center_loglik <- function(p, r){
  # p: |Sc| * 3 matrix, each row is a center
  # r: r = (r1, r2, r3)
  p <- matrix(p, ncol = 3, byrow = FALSE)
  if(nrow(p) == 1){
    loglik = center_loglik(p, r)
  }else{
    loglik <- apply(p, 1, FUN = center_loglik, r = r)
    loglik <- sum(loglik)
  }
  return(loglik)
}

# calculate log acceptance ratio
logacpr <- function(p, r.old, r.new){
  if(is.null(nrow(p))){
    u <- p[1]
    pj <- p[2:4]
  }else{
    u <- p[, 1]
    pj <- p[, 2:4]
  }
  target.old.theta <- sum_center_loglik_theta(u, r.old[1])
  target.new.theta <- sum_center_loglik_theta(u, r.new[1])
  target.old <- sum_center_loglik(pj, r.old) - sum(r.old) + target.old.theta
  target.new <- sum_center_loglik(pj, r.new) - sum(r.new) + target.new.theta
  logjacobian <- log(r.new[1] * r.new[2] * r.new[3] * r.new[4]) - 
    log(r.old[1] * r.old[2] * r.old[3] * r.old[4])
  lacpr <- target.new - target.old + logjacobian
  return(min(0, lacpr))
}

## sampling update for one group
dirichlet.mhsample <- function(p, x, 
                               i, lambda, mu, Sigma){
  d <- length(x)
  # change to log scale
  log.x <- log(x)
  # sample 
  acc <- 0
  iters <- 0
  while(acc == 0){
    iters <- iters + 1
    log.y <- rmvnorm(n = 1, mean = log.x, sigma = lambda * Sigma)
    y <- exp(log.y)
    # acceptance decision
    logalphai <- logacpr(p, x, y)
    u <- runif(1, 0, 1)
    if(log(u) < logalphai){
      x <- y # case of accept
      log.x <- log(x)
      acc <- 1
    }
  }
  
  
  # update
  ri <- 1/sqrt(i+1)
  loglambda <- log(lambda) + ri * (exp(logalphai) - 0.234)
  lambda <- exp(loglambda)
  mu <- mu + ri * (log.x - mu)
  Sigma <- Sigma + 
    ri * (matrix(log.x - mu, ncol = 1) %*% matrix(log.x - mu, ncol = d) - Sigma)
  return(list(x = x, lambda = lambda, mu = mu, Sigma = Sigma, 
              logalphai = logalphai, u = u, y = y, acc = 1/iters))
}
