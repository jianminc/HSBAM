# function to update cluster labels
# calculate dirichlet density
diridense <- function(p, alpha){
  r1 <- alpha[1]
  r2 <- alpha[2]
  r3 <- alpha[3]
  p.pos <- p[1]
  p.neg <- p[2]
  p0 <- p[3]
  loglik <- lgamma(r1 + r2 + r3) - 
    lgamma(r1) - lgamma(r2) - lgamma(r3) + 
    (r1 - 1) * log(p.pos) + (r2 - 1) * log(p.neg) + (r3 - 1) * log(p0) 
  return(exp(loglik))
}

# specify base function here: Gamma(a, a) for all 4 parameters
label.update <- function(nr, j, param.old, uj, pj){
  # get new samples
  new.param.1 <- rgamma(n = nr, shape = 1, rate = 1)
  new.param <- matrix(rgamma(n = nr * 3, 1, 1), nrow = nr)
  new.param <- cbind(new.param.1, new.param)
  # check existing clusters
  param.c <- unique(param.old[-j, 1:4])
  exist <- apply(param.c, 1, FUN = function(x) sum(all.equal(x, param.old[j, 1:4]) == TRUE))
  
  if(sum(exist) == 0){
    new.param[1, ] <- unlist(param.old[j, 1:4])
    new.cluster <- (max(param.old[-j, ]$cluster) + 1):(max(param.old[-j, ]$cluster) + nr - 1)
    new.cluster <- c(param.old$cluster[j], new.cluster)
  }else{
    new.cluster <- (max(param.old[-j, ]$cluster) + 1):(max(param.old[-j, ]$cluster) + nr)
  }
  # sample cluster label
  alpha <- 1 # concentration parameter in DPM
  # allocate clusters
  sc1 <- param.old[-j, ] %>% group_by(cluster) %>% summarise(count = n())
  sc2 <- rep(alpha/nr, nr)
  sc <- c(sc1$count, sc2)
  
  new.param <- cbind(new.param, new.cluster)
  param.all <- rbind(as.matrix(unique(param.old[-j, ])), as.matrix(new.param))
  pu.c <- apply(param.all, 1, FUN = function(x) dgamma(uj, shape = x[1], rate = x[1]))
  #pp.c <- apply(param.all, 1, FUN = function(x) ddirichlet(pj, alpha = x[2:4]))
  pp.c <- apply(param.all, 1, FUN = function(x) diridense(pj, alpha = x[2:4]))
  pi <- pu.c * pp.c * sc
  pi <- pi/(sum(pi))
  ci <- sample(1:nrow(param.all), size = 1, prob = pi, replace = TRUE)
  # use row index to mark clusters: allocate cluster at row ci
  param.j.new <- param.all[ci, ]
  return(param.j.new)
}
