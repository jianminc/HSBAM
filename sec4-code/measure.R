## 
# data.new.pre0: without cap adjustment
# data.new.pre: with cap adjustment

## 1. find exact time that reach the enrollment target
## exact data for both train and test
data.exact <- data1 %>% filter(screen.time < max.time & status %in% c(1, 2)) %>% 
  arrange(screen.time)
data.exact.1 <- data1 %>% filter(screen.time < max.time & status == 1) %>% 
  arrange(screen.time)
data.exact.2 <- data1 %>% filter(screen.time < max.time & status == 2) %>% 
  arrange(screen.time)
date.complete.1 <- data.exact.1$screen.time[target.enroll.sub1]
date.complete.2 <- data.exact.2$screen.time[target.enroll.sub2]
date.complete <- max(date.complete.1, date.complete.2)
exact.date <- c(date.complete.1, date.complete.2, date.complete)

## 2. find prediction for the completion date
##  without adjustment
predate0 <- matrix(0, nrow = n.repeat, ncol = 3)
for (k in 1:n.repeat) {
  data.pre.1 <- data.new.pre0[[k]] %>% filter(status == 1) %>% 
    arrange(screen.time) %>% head(n = target.enroll.sub1)
  predate0[k, 1] <- data.pre.1$screen.time[target.enroll.sub1]
  data.pre.2 <- data.new.pre0[[k]] %>% filter(status == 2) %>% 
    arrange(screen.time) %>% head(n = target.enroll.sub2)
  predate0[k, 2] <- data.pre.2$screen.time[target.enroll.sub2]
  predate0[k, 3] <- max(predate0[k, 1], predate0[k, 2])
}

##  with adjustment
predate <- matrix(0, nrow = n.repeat, ncol = 3)
for (k in 1:n.repeat) {
  data.pre.1 <- data.new.pre[[k]] %>% filter(status == 1) %>% 
    arrange(screen.time) %>% head(n = target.enroll.sub1)
  predate[k, 1] <- data.pre.1$screen.time[target.enroll.sub1]
  data.pre.2 <- data.new.pre[[k]] %>% filter(status == 2) %>% 
    arrange(screen.time) %>% head(n = target.enroll.sub2)
  predate[k, 2] <- data.pre.2$screen.time[target.enroll.sub2]
  predate[k, 3] <- max(predate[k, 1], predate[k, 2])
}

## calculate absolute difference
pre0 <- apply(predate0, 2, median)
pre <- apply(predate, 2, median)
abs.dif1 <- abs(pre0 - exact.date)
abs.dif <- abs(pre - exact.date)

## 3. find prediction within 30, 60, 90 days after the cutoff date

getkdaydif <- function(d){
  # truth
  data.exact <- data1 %>% filter(screen.time <= cutoff + d & status %in% c(1, 2))
  data.exact.1 <- data1 %>% filter(screen.time <= cutoff + d & status == 1)
  data.exact.2 <- data1 %>% filter(screen.time <= cutoff + d & status == 2)
  exact.date <- c(nrow(data.exact.1), nrow(data.exact.2), nrow(data.exact))
  
  # pre without adjust
  predate0 <- matrix(0, nrow = n.repeat, ncol = 3)
  for (k in 1:n.repeat) {
    data.pre.1 <- data.new.pre0[[k]] %>% filter(screen.time <= cutoff + d & status == 1)
    data.pre.2 <- data.new.pre0[[k]] %>% filter(screen.time <= cutoff + d & status == 2)
    predate0[k, 1] <- nrow(data.pre.1)
    predate0[k, 2] <- nrow(data.pre.2)
    predate0[k, 3] <- predate0[k, 1] + predate0[k, 2]
  }
  pre0 <- apply(predate0, 2, median)
  
  # pre with adjust
  predate <- matrix(0, nrow = n.repeat, ncol = 3)
  for (k in 1:n.repeat) {
    data.pre.1 <- data.new.pre[[k]] %>% filter(screen.time <= cutoff + d & status == 1)
    data.pre.2 <- data.new.pre[[k]] %>% filter(screen.time <= cutoff + d & status == 2)
    predate[k, 1] <- nrow(data.pre.1)
    predate[k, 2] <- nrow(data.pre.2)
    predate[k, 3] <- predate[k, 1] + predate[k, 2]
  }
  pre <- apply(predate, 2, median)
  
  abs.dif1 <- abs(pre0 - exact.date)
  abs.dif <- abs(pre - exact.date)
  
  return(list(abs1.dif1 = abs.dif1, abs.dif = abs.dif))
}

# 30 days
res30 <- getkdaydif(30)
res60 <- getkdaydif(60)
res90 <- getkdaydif(90)

