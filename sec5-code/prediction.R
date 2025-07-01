## do prediction
data.new <- vector('list', n.repeat)
data <- data1.final
data <- data %>% filter(period == 0)
nonzero.center <- data %>% filter(censor == 1) %>% 
  group_by(center) %>% summarise(count=n()) %>% 
  filter(count > 0)
data <- data %>% filter(center %in% nonzero.center$center)


# calculate base time before cutoff for candidate centers
# opened center with period 2 IAT
center1 <- center.info %>% filter(plan == 1 & open == 1 & start == 1 & close == 0) %>% 
  dplyr::select(center)
center1.use <- unique(data$center)
center1 <- center1 %>% filter(center %in% center1.use)
group.time <- data1.final %>% filter(screened == 1) %>% 
  group_by(center) %>% summarise(basetime = max(screen.time))
center1 <- center1 %>% left_join(group.time) 
center1$wait.time <- cutoff - center1$basetime
center1$type <- 1
# opened center without period 2 IAT
center2 <- center.info %>% filter(plan == 1 & open == 1 & start == 1 & close == 0) %>% 
  dplyr::select(center)
center2 <- center2 %>% filter(!center %in% center1.use)
if(nrow(center2) > 0){
  center2 <- data1.final %>% filter(center %in% center2$center & censor == 1) %>% 
    mutate(basetime = screen.time, wait.time = cutoff - basetime) %>% 
    dplyr::select(center, basetime, wait.time)
  center2$type <- 2
}
# opened center without patients
center3 <- center.info %>% filter(plan == 1 & open == 1 & start == 0 & close == 0) %>% 
  mutate(basetime = open.time, wait.time = cutoff - basetime) %>% 
  dplyr::select(center, basetime, wait.time)
if(nrow(center3) > 0){
  center3$type <- 3
}
# planned center, not yet opened
max.gap.time.open <- max(data1.final$gap.time.open)
center4 <- center.info %>% filter(plan == 1 & open == 0) %>% 
  mutate(basetime = plan.time, wait.time = cutoff - basetime) %>% 
  dplyr::select(center, basetime, wait.time)
if(nrow(center4) > 0){
  center4$type <- 4
}
# time.open.new <- apply(center3, 1, FUN = function(x){runif(1, x[2], max.time.open)})
#center3 <- data.frame(center = center3$center, basetime = time.open.new, wait.time = 0)
# get all centers
pre.centers <- center1 %>% filter(center %in% center.include)
if(nrow(center2) > 0){
  pre.centers <- rbind(pre.centers, center2)
}
if(nrow(center3) > 0){
  pre.centers <- rbind(pre.centers, center3)
}
if(nrow(center4) > 0){
  pre.centers <- rbind(pre.centers, center4)
}
center.all <- unique(data1$center)
center.train <- unique(data$center)
pre.centers <- pre.centers %>% arrange(center)
# ========================================================================
## sample from posterior
# data.new.enroll <- vector('list', n.repeat)
for(i in 1:n.repeat){
  centers.new <- pre.centers
  new <- matrix(0, nrow = h*J.pre, ncol = 3)
  new[, 3] <- 1
  sampleid <- sample(n.burnin:niter, 1)
  lambda <- iter.data$lambda[sampleid]
  # sample opening time for type 4 centers
  if(sum(pre.centers$type == 4) > 0){
    centers.new$basetime[pre.centers$type == 4] <- 
      apply(center4, 1, FUN = function(x){
        return(x[2] + runif(1, x[3], max(x[3], max.gap.time.open)))
      })
  }
  for(j in seq(J.pre)){
    if(pre.centers$type[j] == 1){ 
      # opened, screened center
      j.id <- which(center.train == pre.centers$center[j])
      uj <- iter.data$u[sampleid, j.id]
      trunc <- c(pre.centers$wait.time[j], rep(0, h - 1))
      pj <- iter.data$p[j.id, , sampleid]
      
      S0 <- exp(-lambda * trunc)
      v <- runif(h, min = 0, max = S0^uj)
      # new IAT
      iat.new <- -log(v) / uj / lambda
      # new screening result
      counts.new <- rmultinom(n = 1, size = h, prob = pj)
      rid <- sample(x = seq(h),  size= h, replace = FALSE)
      status.new <- rep(c(1, 2, 3), time = counts.new)
      status.new <- status.new[rid]
      
      
    }else if(pre.centers$type[j] == 2){
      # opened, only with 1 observation
      j.id <- sample(center.train, size = 1)
      j.id <- which(center.train == j.id)
      thetaj <- iter.data$theta[sampleid, ]
      uj <- rgamma(1, shape = thetaj, rate = thetaj)
      trunc <- c(pre.centers$wait.time[j], rep(0, h-1)) 
      gammaj <- iter.data$gamma[j.id, , sampleid]
      pj <- rdirichlet(1, gammaj)
      tries = 1
      while(sum(is.nan(pj)) > 0 & tries <= 10){
        pj <- rdirichlet(1, gammaj)
        tries <- tries + 1
      }
      if(sum(is.nan(pj)) > 0){
        oneid <- sample(c(1, 2, 3), 1)
        pj <- rep(0, 3)
        pj[oneid] <- 1
      }
      
      
      # new IAT
      S0 <- exp(-lambda * trunc)
      v <- runif(h, min = 0, max = S0^uj)
      # new IAT
      iat.new <- -log(v) / uj / lambda
      # new screening result
      counts.new <- rmultinom(n = 1, size = h, prob = pj)
      rid <- sample(x = seq(h),  size= h, replace = FALSE)
      status.new <- rep(c(1, 2, 3), time = counts.new)
      status.new <- status.new[rid]
      
    }else if(pre.centers$type[j] == 3){ 
      lambda1 <- iter.data$lambda1[sampleid]
      # opened, not yet screened center
      j.id <- sample(center.train, size = 1)
      j.id <- which(center.train == j.id)
      thetaj <- iter.data$theta[sampleid, ]
      uj <- rgamma(1, shape = thetaj, rate = thetaj)
      trunc <- c(pre.centers$wait.time[j], rep(0, h - 1)) 
      gammaj <- iter.data$gamma[j.id, , sampleid]
      pj <- rdirichlet(1, gammaj)
      tries = 1
      while(sum(is.nan(pj)) > 0 & tries <= 10){
        pj <- rdirichlet(1, gammaj)
        tries <- tries + 1
      }
      if(sum(is.nan(pj)) > 0){
        oneid <- sample(c(1, 2, 3), 1)
        pj <- rep(0, 3)
        pj[oneid] <- 1
      }
      
      # new IAT
      S0_1 <- exp(-lambda1 * trunc[1])
      iat.new1 <- -log(runif(1, min = 0, max = S0_1))/lambda1
      S0 <- exp(-lambda * trunc[2:h])
      v <- runif(h-1, min = 0, max = S0^uj)
      iat.new <- c(iat.new1, -log(v) / uj / lambda)
      # new screening result
      counts.new <- rmultinom(n = 1, size = h, prob = pj)
      rid <- sample(x = seq(h),  size= h, replace = FALSE)
      status.new <- rep(c(1, 2, 3), time = counts.new)
      status.new <- status.new[rid]
      
      
    }else{ 
      # planned, not opened yet center
      lambda1 <- iter.data$lambda1[sampleid]
      j.id <- sample(center.train, size = 1)
      j.id <- which(center.train == j.id)
      # theta.sampled <- iter.data$theta[sampleid, ]
      thetaj <- iter.data$theta[sampleid, j.id]
      #thetaj <- sample(theta.sampled, size = 1)
      uj <- rgamma(1, shape = thetaj, rate = thetaj)
      trunc <- rep(0, h) 
      gammaj <- iter.data$gamma[j.id, , sampleid]
      pj <- rdirichlet(1, gammaj)
      tries = 1
      while(sum(is.nan(pj)) > 0 & tries <= 10){
        pj <- rdirichlet(1, gammaj)
        tries <- tries + 1
      }
      if(sum(is.nan(pj)) > 0){
        oneid <- sample(c(1, 2, 3), 1)
        pj <- rep(0, 3)
        pj[oneid] <- 1
      }
      
      # new IAT
      S0_1 <- exp(-lambda1 * trunc[1])
      iat.new1 <- -log(runif(1, min = 0, max = S0_1))/lambda1
      S0 <- exp(-lambda * trunc[2:h])
      v <- runif(h-1, min = 0, max = S0^uj)
      iat.new <-c(iat.new1, -log(v) / uj / lambda)
      # new screening result
      counts.new <- rmultinom(n = 1, size = h, prob = pj)
      rid <- sample(x = seq(h),  size= h, replace = FALSE)
      status.new <- rep(c(1, 2, 3), time = counts.new)
      status.new <- status.new[rid]
    }
    new[((j-1)*h+1):(j*h), 1] <- iat.new
    new[((j-1)*h+1):(j*h), 2] <- status.new
    
  }
  new <- as.data.frame(new)
  colnames(new) <- c('iat.new', 'status', 'screened')
  new$center <- rep(centers.new$center, each = h)
  time.exact.new <- new %>% group_by(center) %>% 
    summarise(time = cumsum(iat.new))
  new$screen.time <- time.exact.new$time
  new <- new %>% left_join(centers.new, by = c('center'))
  new$screen.time <- new$screen.time + new$basetime
  if(sum(new$screen.time < cutoff) > 0){
    break
  }
  new$enroll.window <- 0
  new.n.enroll <- sum(new$status %in% c(1, 2))
  #new$enroll.window[new$status %in% c(1, 2)] <- runif(n = new.n.enroll, 1, 7)
  new$enroll.time <- new$screen.time + new$enroll.window
  new$source <- 'predict'
  new <- new %>% arrange(screen.time) 
  new <- new %>% dplyr::select(status, screened, 
                               screen.time, enroll.time, source,center)
  ## adjustment for screened, not enrolled cases
  if(length(preenroll.id) > 0){
    enroll.wait.time <- cutoff - data1.train.preenroll$screen.time
    gap.time.enroll <- runif(n = length(preenroll.id), 
                             min = enroll.wait.time, 
                             max = 7)
    preenroll.enroll.time <- data1.train.preenroll$screen.time + gap.time.enroll
  }else{
    preenroll.enroll.time <- NULL
  }
  data.new[[i]] <- list(new = new, preenroll.enroll.time = preenroll.enroll.time)
}



## get prediction (also for plots)
plotdata <- data1.train %>% 
  dplyr::select(status, screened, screen.time, 
                enroll.time, center) %>% mutate(source = 'truth')
## ==========================================================
data.new.pre <- vector('list', n.repeat)
plotdata.new <- plotdata
if(length(preenroll.id) > 0){
  for (k in 1:n.repeat) {
    plotdata.new[preenroll.id, ]$enroll.time <- data.new[[k]]$preenroll.enroll.time
    data.new.pre[[k]] <- data.new[[k]]$new %>% filter(screen.time <= cutoff + max.time * 10)
    data.new.pre[[k]] <- rbind(plotdata.new, data.new.pre[[k]])
  }
}else{
  for (k in 1:n.repeat) {
    data.new.pre[[k]] <- data.new[[k]]$new %>% filter(screen.time <= cutoff + max.time * 10)
    data.new.pre[[k]] <- rbind(plotdata, data.new.pre[[k]])
  }
}


# adjust the predictions by target enrollment for each disease sub-type in each region
# stop enrollment when the cap is reached
data.new.pre0 <- data.new.pre
data.new.pre.adj <- lapply(data.new.pre, FUN = function(x){
  finish.time1.type1 <- x %>% filter(center %in% center.list1 & status == 1) %>% 
    arrange(enroll.time) %>% slice(target1.1.min:target1.1.min) %>% 
    dplyr::select(enroll.time)
  finish.time1.type2 <- x %>% filter(center %in% center.list1 & status == 2) %>% 
    arrange(enroll.time) %>% slice(target1.2.min:target1.2.min) %>% 
    dplyr::select(enroll.time)
  if(nrow(finish.time1.type1) > 0 & nrow(finish.time1.type2) > 0){
    finish.time1 <- max(finish.time1.type1$enroll.time, finish.time1.type2$enroll.time)
  }else{
    finish.time1 <- max.time
    #print('sample number insufficient')
    #break
  }
  
  finish.time2.type1 <- x %>% filter(center %in% center.list2 & status == 1) %>% 
    arrange(enroll.time) %>% slice(target2.1.min:target2.1.min) %>% 
    dplyr::select(enroll.time)
  finish.time2.type2 <- x %>% filter(center %in% center.list2 & status == 2) %>% 
    arrange(enroll.time) %>% slice(target2.2.min:target2.2.min) %>% 
    dplyr::select(enroll.time)
  if(nrow(finish.time2.type1) > 0 & nrow(finish.time2.type2) > 0){
    finish.time2 <- max(finish.time2.type1$enroll.time, finish.time2.type2$enroll.time)
  }else{
    finish.time2 <- max.time
    #print('sample number insufficient')
    #break
  }
  
  finish.time3.type1 <- x %>% filter(center %in% center.list3 & status == 1) %>% 
    arrange(enroll.time) %>% slice(target3.1.min:target3.1.min) %>% 
    dplyr::select(enroll.time)
  finish.time3.type2 <- x %>% filter(center %in% center.list3 & status == 2) %>% 
    arrange(enroll.time) %>% slice(target3.2.min:target3.2.min) %>% 
    dplyr::select(enroll.time)
  if(nrow(finish.time3.type1) > 0 & nrow(finish.time3.type2) > 0){
    finish.time3 <- max(finish.time3.type1$enroll.time, finish.time3.type2$enroll.time)
  }else{
    finish.time3 <- max.time
    #print('sample number insufficient')
    #break
  }
  
  x <- x %>% filter(!(center %in% center.list1 & screen.time >= finish.time1)) 
  x <- x %>% filter(!(center %in% center.list2 & screen.time >= finish.time2)) 
  x <- x %>% filter(!(center %in% center.list3 & screen.time >= finish.time3)) 
  
  # region 1
  n1.1 <- nrow(x %>% filter(center %in% center.list1 & status == 1))
  if(n1.1 > target1.1.max){
    x <- rmdata(x, center.list1, target1.1.max, 1)
  }
  n1.2 <- nrow(x %>% filter(center %in% center.list1 & status == 2))
  if(n1.2 > target1.2.max){
    x <- rmdata(x, center.list1, target1.2.max, 2)
  }
  
  # region 2
  n2.1 <- nrow(x %>% filter(center %in% center.list2 & status == 1))
  if(n2.1 > target2.1.max){
    x <- rmdata(x, center.list2, target2.1.max, 1)
  }
  n2.2 <- nrow(x %>% filter(center %in% center.list2 & status == 2))
  if(n2.2 > target2.2.max){
    x <- rmdata(x, center.list2, target2.2.max, 2)
  }
  
  # region 3
  n3.1 <- nrow(x %>% filter(center %in% center.list3 & status == 1))
  if(n3.1 > target3.1.max){
    x <- rmdata(x, center.list3, target3.1.max, 1)
  }
  n3.2 <- nrow(x %>% filter(center %in% center.list3 & status == 2))
  if(n3.2 > target3.2.max){
    x <- rmdata(x, center.list3, target3.2.max, 2)
  }
  
  return(x)
})

data.new.pre <- data.new.pre.adj

# get estimated region enrollment completion time
time <- lapply(data.new.pre.adj, FUN = function(x){
  finish.time1.type1 <- x %>% filter(center %in% center.list1 & status == 1) %>% 
    arrange(enroll.time) %>% slice(target1.1.min:target1.1.min) %>% 
    dplyr::select(enroll.time)
  finish.time1.type2 <- x %>% filter(center %in% center.list1 & status == 2) %>% 
    arrange(enroll.time) %>% slice(target1.2.min:target1.2.min) %>% 
    dplyr::select(enroll.time)
  if(nrow(finish.time1.type1) > 0 & nrow(finish.time1.type2) > 0){
    finish.time1 <- max(finish.time1.type1$enroll.time, finish.time1.type2$enroll.time)
  }else{
    finish.time1 <- max.time
    #print('sample number insufficient')
    #break
  }
  
  finish.time2.type1 <- x %>% filter(center %in% center.list2 & status == 1) %>% 
    arrange(enroll.time) %>% slice(target2.1.min:target2.1.min) %>% 
    dplyr::select(enroll.time)
  finish.time2.type2 <- x %>% filter(center %in% center.list2 & status == 2) %>% 
    arrange(enroll.time) %>% slice(target2.2.min:target2.2.min) %>% 
    dplyr::select(enroll.time)
  if(nrow(finish.time2.type1) > 0 & nrow(finish.time2.type2) > 0){
    finish.time2 <- max(finish.time2.type1$enroll.time, finish.time2.type2$enroll.time)
  }else{
    finish.time2 <- max.time
    #print('sample number insufficient')
    #break
  }
  
  finish.time3.type1 <- x %>% filter(center %in% center.list3 & status == 1) %>% 
    arrange(enroll.time) %>% slice(target3.1.min:target3.1.min) %>% 
    dplyr::select(enroll.time)
  finish.time3.type2 <- x %>% filter(center %in% center.list3 & status == 2) %>% 
    arrange(enroll.time) %>% slice(target3.2.min:target3.2.min) %>% 
    dplyr::select(enroll.time)
  if(nrow(finish.time3.type1) > 0 & nrow(finish.time3.type2) > 0){
    finish.time3 <- max(finish.time3.type1$enroll.time, finish.time3.type2$enroll.time)
  }else{
    finish.time3 <- max.time
    #print('sample number insufficient')
    #break
  }
  
  finish.time4.type1 <- x %>% filter(status == 1) %>% 
    arrange(enroll.time) %>% slice(target.enroll.sub1:target.enroll.sub1) %>% 
    dplyr::select(enroll.time)
  finish.time4.type2 <- x %>% filter(status == 2) %>% 
    arrange(enroll.time) %>% slice(target.enroll.sub2:target.enroll.sub2) %>% 
    dplyr::select(enroll.time)
  finish.time4 <- max(finish.time4.type1, finish.time4.type2, 
                      finish.time1, finish.time2, finish.time3)
  
  return(c(finish.time1, finish.time2, finish.time3, finish.time4))
})

t1 <- sapply(time, function(x) x[1])
c(median(t1), quantile(t1, c(0.025,0.975))) / 30
t2 <- sapply(time, function(x) x[2])
c(median(t2), quantile(t2, c(0.025,0.975))) / 30
t3 <- sapply(time, function(x) x[3])
c(median(t3), quantile(t3, c(0.025,0.975))) / 30
t4 <- sapply(time, function(x) x[4])
c(median(t4), quantile(t4, c(0.025,0.975))) / 30


# true completion time for each region
max(data1 %>% filter(status == 1 & center %in% center.list1) %>% arrange(enroll.time) %>% 
    dplyr::select(enroll.time) %>% slice(target1.1.min:target1.1.min) / 30, 
    data1 %>% filter(status == 2 & center %in% center.list1) %>% arrange(enroll.time) %>% 
    dplyr::select(enroll.time) %>% slice(target1.2.min:target1.2.min) / 30)

max(data1 %>% filter(status == 1 & center %in% center.list2) %>% arrange(enroll.time) %>% 
      dplyr::select(enroll.time) %>% slice(target2.1.min:target2.1.min) / 30, 
    data1 %>% filter(status == 2 & center %in% center.list2) %>% arrange(enroll.time) %>% 
      dplyr::select(enroll.time) %>% slice(target2.2.min:target2.2.min) / 30)

max(data1 %>% filter(status == 1 & center %in% center.list3) %>% arrange(enroll.time) %>% 
      dplyr::select(enroll.time) %>% slice(target3.1.min:target3.1.min) / 30, 
    data1 %>% filter(status == 2 & center %in% center.list3) %>% arrange(enroll.time) %>% 
      dplyr::select(enroll.time) %>% slice(target3.2.min:target3.2.min) / 30)

max(data1 %>% filter(status == 1) %>% arrange(enroll.time) %>% 
      dplyr::select(enroll.time) %>% slice(target.enroll.sub1:target.enroll.sub1) / 30, 
    data1 %>% filter(status == 2) %>% arrange(enroll.time) %>% 
      dplyr::select(enroll.time) %>% slice(target.enroll.sub2:target.enroll.sub2) / 30)

