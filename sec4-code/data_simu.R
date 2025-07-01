rmdata <- function(data, centerlist, target, st){
  t.max <-  data %>% filter(center %in% centerlist & status == st) %>% 
    arrange(enroll.time) %>% slice(target:target) %>% dplyr::select(enroll.time)
  data <- data %>% filter(!(center %in% centerlist & 
                              status == st & 
                              screen.time > t.max$enroll.time))
  return(data)
}


# data_simu <- function(K = 3, nc = c(25, 15, 10), npt = 50,
#                       thetak = c(0.5, 1, 6), rk= c(2, 1, 0.5, 2, 4, 1, 1, 1, 1),
#                       lambda1 = 0.05, lambda2 = 0.1, 
#                       cutoff = 300, max.time = 365 * 3){
#   
  
  # total number of centers
  J <- sum(nc)
  n <- rep(npt, J)
  # stage 0
  #=======================================================
  # sample planning date for centers: assume 5 centers are planned at time 0
  # assume the planned date for each center is known at the beginning of the study
  plan.time <- floor(runif(J, 0, 180))
  plan.time[sample(1:J, 5, replace = FALSE)] <- c(rep(0, 5))
  
  # stage 1
  #=======================================================
  # sample time from planning to opening
  gap.time.open <- floor(runif(J, 0, 60))
  
  # stage 2
  #=======================================================
  # sample first IAT at each center
  # generate iat for the first enrollment
  iat1 <- rgamma(n = J, shape = 1, 
                 rate = lambda1)
  iat1 <- ceiling(iat1)
  data0 <- data.frame(center = 1:J, iat = iat1, stage = 1)
  # sample first outcome
  
  
  
  # stage 3
  #========================================================
  # consecutive inter-arrival times
  #=========================================================
  # cluster specific parameter for mixture gamma and mixture dirichlet distribution
  rk <- matrix(rk, byrow = TRUE, nrow = 3)
  #=========================================================
  # center specific parameter
  params.center <- matrix(0, nrow = sum(nc), ncol = 4)
  id <- c(0, cumsum(nc))
  for(i in seq(K)){
    uj <- rgamma(n = nc[i], shape = thetak[i], rate = thetak[i])
    #uj <- rgamma(n = nc[i], shape = thetak[i], rate = 1)
    pj <- rdirichlet(n = nc[i], alpha = rk[i, ])
    param <- cbind(uj, pj)
    params.center[(id[i] + 1):(id[i + 1]), ] <- param
  }
  #========================================================
  # IAT in each center
  id <- c(0, cumsum(n))
  data1 <- matrix(0, nrow = sum(n), ncol = 1)
  # generate IAT for each patient
  for(i in seq(J)){
    tij <- rgamma(n = n[i], shape = 1, 
                  rate = lambda2 * params.center[i, 1])
    data1[(id[i] + 1):(id[i + 1]), 1] <- tij
  }
  data1 <- data.frame(time = data1[, 1])
  data1$center <- rep(1:J, times = n)
  
  # exclude extreme IAT values
  data1$exclude <- as.numeric(data1$time > quantile(data1$time, 0.8))
  n.adj <- data1 %>% group_by(center) %>% summarise(n.exc = sum(exclude))
  n <- n - n.adj$n.exc
  data1 <- data1 %>% filter(exclude == 0) %>% dplyr::select(-exclude)
  # convert time into integers: convert to date scale
  data1$iat <- floor(data1$time) + 1
  data1 <- data1 %>% dplyr::select(-time)
  data1$stage <- 2
  
  # combine iat1 and the other iats
  data1 <- rbind(data0, data1) %>% arrange(center, stage)
  n <- n + 1
  # # drop empty centers after adjustment
  # if(sum(n == 0) > 0){
  #   n <- n[n > 0]
  #   params.center <- params.center[n > 0, ]
  #   J <- length(n)
  # }
  
  # get cumulative time: screened time
  time.exact <- data1 %>% group_by(center) %>% summarise(time = cumsum(iat))
  data1$screen.time <- time.exact$time
  
  #=================================================================
  # get center time gap before close time
  center.close.info <- matrix(0, nrow = J, ncol = 3)
  centers <- unique(data1$center)
  for(i in 1:J){
    if(n[i] > 1){
      closei <- data1 %>% filter(center == centers[i]) %>% 
        group_by(center) %>% 
        arrange(desc(screen.time)) %>% 
        summarise(mtime = head(screen.time, n = 2)[1] - head(screen.time, n = 2)[2], 
                  gap.time.close = round(runif(1, 0, mtime)) + 1) %>% 
        dplyr::select(center, mtime, gap.time.close)
      center.close.info[i, ] <- unlist(closei)
    }else{
      mtime <- runif(1, 90, 360)
      center.close.info[i, ] <- c(centers[i], mtime, mtime)
    }
  }
  center.close.info <- as.data.frame(center.close.info)
  names(center.close.info) <- c('center', 'mtime', 'gap.time.close')
  
  #=======================================================
  # get all the center level information
  center.info <- data.frame(center = unique(data1$center), 
                            gap.time.open = gap.time.open, 
                            gap.time.close = center.close.info$gap.time.close, 
                            plan.time = plan.time, 
                            open.time = plan.time + gap.time.open)
  center.info$open.time <- center.info$plan.time + center.info$gap.time.open
  data1 <- data1 %>% left_join(center.info, by = c('center'))
  data1$screen.time <- data1$screen.time + data1$gap.time.open + data1$plan.time
  close.time <- data1 %>% group_by(center) %>% 
    summarise(close.time = max(screen.time)) 
  center.info$close.time <- close.time$close.time + center.info$gap.time.close
  
  #=============================================
  # define already planned centers
  center.info$plan <- as.numeric(center.info$plan.time < cutoff) 
  # define already opened centers
  center.info$open <- as.numeric(center.info$open.time < cutoff)
  # define closed centers 
  center.info$close <- as.numeric(center.info$close.time < cutoff)
  # define already started screening centers
  data1$screened <- as.numeric(data1$screen.time < cutoff)
  center.start <- data1 %>% filter(screened == 1) %>% 
    group_by(center) %>% summarise(count = n())
  center.start <- data.frame(center = unique(data1$center), start = as.numeric(unique(data1$center) %in% center.start$center))
  center.info$start <- center.start$start 
  
  # ==================================================================================
  # get screening result
  # status: 1 subtype 1, enrolled
  # status: 2 subtype 2, enrolled
  # status: 3 not enrolled/screening failure case
  data1$status <- 0
  id <- c(0, cumsum(n))
  # get biomarker testing result
  for(i in seq(J)){
    counts <- rmultinom(n = 1, size = n[i], prob = params.center[i, 2:4])
    status <- rep(c(1, 2, 3), times = counts)
    # randomly assign status
    rid <- sample(seq(n[i]), n[i], replace = FALSE)
    data1$status[(id[i] + 1):(id[i + 1])] <- status[rid]
  }
  
  #=================================================================
  # add enrollment time window: in integer
  n.enroll <- sum(data1$status %in% c(1, 2))
  data1$gap.time.enroll <- 0
  #data1$gap.time.enroll[data1$status %in% c(1, 2)] <- round(runif(n.enroll, 1, 7))
  data1$enroll.time <- data1$screen.time + data1$gap.time.enroll
  data1$enroll.time[data1$status == 3] <- 0
  
  #==================================================================
  # adjust closing time and data according to pre specified target enrollment info
  # first 1-10 centers: recruit up to 50 pts (in the first cluster)
  # second 16-20 centers: recruit up to 50 pts (in the second cluster)
   center.list1 <- c(1:10)[1:10 %in% center.info$center]
   center.list2 <- c(11:25)[11:25 %in% center.info$center]
   center.list3 <- c(26:40)[26:40 %in% center.info$center]
   center.list4 <- c(41:50)[41:50 %in% center.info$center]
  
  #center.list1 <- c(1:7)[1:8 %in% center.info$center]
  #center.list2 <- c(8:15)[8:15 %in% center.info$center]
  #center.list3 <- c(16:25)[16:25 %in% center.info$center]
  #center.list4 <- c(26:30)[26:30 %in% center.info$center]
  # 
  # get exact finish time for each center, based on minimum requirements
  finish.time1.type1 <- data1 %>% filter(center %in% center.list1 & status == 1) %>% 
    arrange(enroll.time) %>% slice(target1.1.min:target1.1.min) %>% 
    dplyr::select(enroll.time)
  finish.time1.type2 <- data1 %>% filter(center %in% center.list1 & status == 2) %>% 
    arrange(enroll.time) %>% slice(target1.2.min:target1.2.min) %>% 
    dplyr::select(enroll.time)
  if(nrow(finish.time1.type1) > 0 & nrow(finish.time1.type2) > 0){
    finish.time1 <- max(finish.time1.type1$enroll.time, finish.time1.type2$enroll.time)
  }else{
    #print('sample number insufficient')
    #break
    finish.time1 <- max.time
  }
  
  finish.time2.type1 <- data1 %>% filter(center %in% center.list2 & status == 1) %>% 
    arrange(enroll.time) %>% slice(target2.1.min:target2.1.min) %>% 
    dplyr::select(enroll.time)
  finish.time2.type2 <- data1 %>% filter(center %in% center.list2 & status == 2) %>% 
    arrange(enroll.time) %>% slice(target2.2.min:target2.2.min) %>% 
    dplyr::select(enroll.time)
  if(nrow(finish.time2.type1) > 0 & nrow(finish.time2.type2) > 0){
    finish.time2 <- max(finish.time2.type1$enroll.time, finish.time2.type2$enroll.time)
  }else{
    #print('sample number insufficient')
    #break
    finish.time2 <- max.time
  }
  
  finish.time3.type1 <- data1 %>% filter(center %in% center.list3 & status == 1) %>% 
    arrange(enroll.time) %>% slice(target3.1.min:target3.1.min) %>% 
    dplyr::select(enroll.time)
  finish.time3.type2 <- data1 %>% filter(center %in% center.list3 & status == 2) %>% 
    arrange(enroll.time) %>% slice(target3.2.min:target3.2.min) %>% 
    dplyr::select(enroll.time)
  if(nrow(finish.time3.type1) > 0 & nrow(finish.time3.type2) > 0){
    finish.time3 <- max(finish.time3.type1$enroll.time, finish.time3.type2$enroll.time)
  }else{
    #print('sample number insufficient')
    #break
    finish.time3 <- max.time
  }
  
  # adjust center closing time
  # center.info$close.time[center.info$center %in% center.list1] <- 
  #   pmin(center.info$close.time[center.info$center %in% center.list1] ,finish.time1)
  # center.info$close.time[center.info$center %in% center.list2] <- 
  #   pmin(center.info$close.time[center.info$center %in% center.list2] ,finish.time2)
  # center.info$close.time[center.info$center %in% center.list3] <- 
  #   pmin(center.info$close.time[center.info$center %in% center.list3] ,finish.time3)
  # #center.info$close.time[center.info$center %in% center.list4] <- 
  #  pmin(center.info$close.time[center.info$center %in% center.list4] ,finish.time4)
  center.info$close <- as.numeric(center.info$close.time < cutoff)
  
  
  # remove observations from the whole data set after closing date
  # data1 <- data1 %>% filter(!(center %in% center.list1 & screen.time >= finish.time1)) 
  # data1 <- data1 %>% filter(!(center %in% center.list2 & screen.time >= finish.time2)) 
  # data1 <- data1 %>% filter(!(center %in% center.list3 & screen.time >= finish.time3)) 
  # #data1 <- data1 %>% filter(!(center %in% center.list4 & screen.time >= finish.time4)) 
  
  # do adjustments based on target.max
  # region 1
  n1.1 <- nrow(data1 %>% filter(center %in% center.list1 & status == 1))
  if(n1.1 > target1.1.max){
    data1 <- rmdata(data1, center.list1, target1.1.max, 1)
  }
  n1.2 <- nrow(data1 %>% filter(center %in% center.list1 & status == 2))
  if(n1.2 > target1.2.max){
    data1 <- rmdata(data1, center.list1, target1.2.max, 2)
  }
  
  # region 2
  n2.1 <- nrow(data1 %>% filter(center %in% center.list2 & status == 1))
  if(n2.1 > target2.1.max){
    data1 <- rmdata(data1, center.list2, target2.1.max, 1)
  }
  n2.2 <- nrow(data1 %>% filter(center %in% center.list2 & status == 2))
  if(n2.2 > target2.2.max){
    data1 <- rmdata(data1, center.list2, target2.2.max, 2)
  }
  
  # region 3
  n3.1 <- nrow(data1 %>% filter(center %in% center.list3 & status == 1))
  if(n3.1 > target3.1.max){
    data1 <- rmdata(data1, center.list3, target3.1.max, 1)
  }
  n3.2 <- nrow(data1 %>% filter(center %in% center.list3 & status == 2))
  if(n3.2 > target3.2.max){
    data1 <- rmdata(data1, center.list3, target3.2.max, 2)
  }
  
  #==================================================================
  # adjust center info after adjustment
  center.start <- data1 %>% filter(screened == 1) %>% 
    group_by(center) %>% summarise(count = n())
  center.start <- setdiff(center.info$center, center.start$center)
  center.info$start[center.info$center %in% center.start] <- 0
 
  
  cap.center <- data1 %>% filter(screen.time >= cutoff) %>% group_by(center) %>% 
    summarise(count = n())
  cap.center <- cap.center$center[cap.center$count == 0]
  if(length(cap.center > 0)){
    center.info$close[center.info$center %in% cap.center] <- 1
  }
  
  #==================================================================
  # get train and test data
  data1 <- data1 %>% dplyr::select(screened, center, 
                                   iat, gap.time.open, gap.time.enroll, 
                                   plan.time, open.time, 
                                   screen.time, enroll.time, status) %>% 
    mutate(censor = 0)
  data1$enroll <- as.numeric(data1$enroll.time < cutoff & data1$enroll.time > 0)
  # get IAT period label
  period <- data1 %>% group_by(center) %>% 
    mutate(period = as.numeric(row_number() == 1))
  data1$period <- period$period
  
  data1.train <- data1 %>% filter(screened == 1)
  # exclude tail observations to mimic real situation: consider a period of 3 year
  data1.test <- data1 %>% filter(screened == 0 & screen.time < max.time)
  
  #===================================================================
  # add waiting time information to the training data
  data1.train$censor <- 1
  data1.test$censor <- 1
  # calculate waiting time for screening process
  # center: already start screening 
  screen.time.wait <- data1 %>% filter(screened == 1) %>% 
    group_by(center) %>% 
    summarise(time.wait = cutoff - max(screen.time), 
              period = 0)
  # adjust waiting time for closed center
  screen.time.wait1 <- center.info %>% filter(close == 1) %>% 
    dplyr::select(center, gap.time.close)
  if(nrow(screen.time.wait1) > 0){
    screen.time.wait <- screen.time.wait %>% 
      left_join(screen.time.wait1, by = c('center'))
    adj.id <- which(!is.na(screen.time.wait$gap.time.close))
    screen.time.wait$time.wait[adj.id] <- screen.time.wait$gap.time.close[adj.id]
    screen.time.wait <- screen.time.wait %>% dplyr::select(-gap.time.close)
  }
  # opened not yet screened center
  center.id <- center.info$center[center.info$open == 1 & center.info$start == 0]
  if(length(center.id) > 0){
    screen.time.wait2 <- data1[!duplicated(data1$center), ] %>% 
      filter(center %in% center.id) %>% 
      mutate(time.wait = cutoff - open.time, period = 1) %>%
      dplyr::select(center, time.wait, period)
    screen.time.wait <- rbind(screen.time.wait, screen.time.wait2)
  }
  screen.time.wait <- screen.time.wait %>% left_join(center.info, by = c('center'))
  screen.time.wait <- data.frame(screened = 0,
                                 center = screen.time.wait$center, 
                                 iat = screen.time.wait$time.wait,
                                 gap.time.open = screen.time.wait$gap.time.open, 
                                 gap.time.enroll = 0, 
                                 plan.time = screen.time.wait$plan.time, 
                                 open.time = screen.time.wait$open.time,
                                 screen.time = 0, 
                                 enroll.time = 0,
                                 status = 0, 
                                 censor = 0, 
                                 enroll = 0, 
                                 period = screen.time.wait$period)
  data1.final <- rbind(data1.train, screen.time.wait)
  #===========================================================================
  ## waiting time for enrollment
  ## get screened but not yet enrolled observations
  data1.train.preenroll <- data1.train %>% 
    filter(status %in% c(1, 2) & enroll == 0)
  preenroll.id <- which(data1.train$status %in% c(1, 2) & data1.train$enroll == 0)
  
  
  
  
#   return(list())
#   
# }