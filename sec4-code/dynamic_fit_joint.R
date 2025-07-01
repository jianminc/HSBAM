# folded in a function
jointfit <- function(data1.final, niter){
  # set initial values
  data <- data1.final
  J0 <- length(unique(data$center))
  
  ## data for period 1
  n.censor.t1 <- nrow(data %>% filter(period == 1 & censor == 1))
  time.1j <- sum(data$iat[data$period == 1])
  
  ## counts for period 1 and 2 combined
  n.pos <- data %>% group_by(center) %>% summarise(ncenter = sum(status == 1))
  n.neg <- data %>% group_by(center) %>% summarise(ncenter = sum(status == 2))
  n.failed <- data %>% group_by(center) %>% summarise(ncenter = sum(status == 3))
  n.j <- data.frame(center = n.pos$center, pos = n.pos$ncenter, 
                    neg = n.neg$ncenter, failed = n.failed$ncenter)
  
  
  ## data for period 2
  data <- data %>% filter(period == 0)
  # only use center having IAT in period2
  nonzero.center <- data %>% filter(censor == 1) %>% 
    group_by(center) %>% summarise(count=n()) %>% 
    filter(count > 0)
  data <- data %>% filter(center %in% nonzero.center$center)
  n.j <- n.j %>% filter(center %in% nonzero.center$center)
  n.j <- as.matrix(n.j[, 2:4])
  
  
  J <- length(unique(data$center))
  u0 <- sample(1:10, 1)
  u <- rep(u0, J)
  p0 <- rdirichlet(1, c(1, 1, 1))
  p <- matrix(p0, nrow = J, ncol = 3)
  lambda <- sample(1:10, 1)
  theta <- rep(sample(1:10, 1), J)
  gamma <- matrix(sample(1:10, 1), nrow = J, ncol = 3)
  
  # calculate intermediate values from data to avoid repeated computation
  n.censor <- data %>% group_by(center) %>% summarise(ncensor = sum(censor == 1))
  n.censor.all <- sum(n.censor$ncensor)
  time.j <- data %>% group_by(center) %>% summarise(time = sum(iat))
  
  # remove center
  iter.data <- list(lambda1 = rep(0, niter),
                    lambda = matrix(0, nrow = niter), 
                    u = matrix(0, nrow = niter, ncol = J), 
                    p = array(0, dim = c(J, 3, niter)), 
                    theta = matrix(0, nrow = niter, ncol = J),
                    gamma = array(0, dim = c(J, 3, niter))
  )
  
  # set initial cluster parameters and cluster label
  data.center <- matrix(0, nrow = J, ncol = 5)
  # 1st col: theta
  data.center[, 1] <- theta
  # 2rd-4th col: gamma
  data.center[, 2:4] <- gamma
  # 5th col: label
  data.center[, 5] <- 1:J
  data.center <- data.frame(data.center)
  colnames(data.center)[1:5] <- c('theta', 'gamma1', 'gamma2', 'gamma3', 'cluster')
  
  # also save AMH parameters in data.amh matrix
  nu2 <- rep(2.38^2/4, J) # lambda in AMH
  # for r and theta together
  mu <- matrix(0, nrow = J, ncol = 4)
  Sigma <- as.vector(diag(c(1, 1, 1, 1)))
  Sigma <- matrix(Sigma, nrow = J, ncol = 16, byrow = TRUE)
  data.amh <- cbind(nu2, mu, Sigma)
  data.amh.0 <- data.amh[1, ]
  acc <- 0
  for(i in seq(niter)){
    # sample lambda 1
    iter.data$lambda1[i] <- lambda1_Gibbs(n.censor.t1, time.1j) 
    # update lambda
    lambda <- lambda_Gibbs(n.censor.all, time.j, u)
    iter.data$lambda[i] <- lambda
    # update u for each center
    u <- u_Gibbs(n.censor, time.j, lambda, theta)
    iter.data$u[i, ] <- u
    # update p for each center
    p <- p_Gibbs(n.j, gamma)
    if(sum(p < 1e-20) > 0){
      p[which(p < 1e-20, arr.ind = TRUE)] <- 10^-20
      sump <- apply(p, 1, sum)
      p <- diag(1/sump) %*% p 
    }
    iter.data$p[, , i] <- p
    # update cluster label
    nr <- 3 # >= 2
    for(j in seq(J)){
      param.j.new <- label.update(nr, j, data.center, u[j], p[j, ])
      # amh parameter need to follow the change of cluster index
      exist.cluster <- which(data.center$cluster == param.j.new[5])
      if(length(exist.cluster) == 1){
        data.amh[j, ] <- data.amh[data.center$cluster == param.j.new[5], ]
      }else if(length(exist.cluster) > 1){
        data.amh[j, ] <- unique(data.amh[data.center$cluster == param.j.new[5], ])
      }else{
        data.amh[j, ] <- data.amh.0
      }
      data.center[j, ] <- param.j.new
      
    }
    # update cluster parameters: only non-empty cluster will be updated
    for(s in unique(data.center$cluster)){
      data.center.s <- unique(data.center %>% filter(cluster == s))
      nsc <- sum(data.center$cluster == s)
      if(nsc == 1){
        data.amh.s <- data.amh[data.center$cluster == s, ]
        p.s <- c(u[which(data.center$cluster == s)], 
                 p[which(data.center$cluster == s), ])
      }else{
        data.amh.s <- unique(data.amh[data.center$cluster == s, ])
        u.s <- u[which(data.center$cluster == s)]
        p.s <- p[which(data.center$cluster == s), ]
        p.s <- cbind(u.s, p.s)
      }
      
      # for r and theta together
      r.s <- unlist(data.center.s)[1:4]
      nu <- data.amh.s[1]       # lambda parameter in AMH
      mu <- data.amh.s[2:5]     # prior mean
      Sigma <- matrix(data.amh.s[6:21], 4, 4)     # prior covariance matrix 
      samples.r <- dirichlet.mhsample(p = p.s, x = r.s, i = i, 
                                      lambda = nu, mu = mu, Sigma = Sigma)
      
      
      # update theta and r
      theta.new <- samples.r$y[1]
      r.new <- samples.r$y[2:4]
      data.center[which(data.center$cluster == s), 1:4] <- matrix(c(theta.new, r.new), nrow = nsc, ncol = 4, byrow = TRUE)
      data.amh[which(data.center$cluster == s), ] <- matrix(c(samples.r$lambda, samples.r$mu, as.vector(samples.r$Sigma)), 
                                                            nrow = nsc, ncol = 21, byrow = TRUE)
      acc <- c(acc, samples.r$acc)
    }
    iter.data$theta[i, ] <- as.matrix(data.center[, 1])
    iter.data$gamma[ , , i] <- as.matrix(data.center[, 2:4])
    theta <- as.matrix(data.center[, 1])
    gamma <- as.matrix(data.center[, 2:4])
    #print(i)
  }
  return(iter.data)
}
