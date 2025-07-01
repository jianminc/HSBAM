## Code to conduct simulation studies in Section 4. This is the code for just one 
## repetition. 

source('preamble.R')

## =============================================================================
## Step 1 : generate simulated data 
seed <- 1
set.seed(seed)
# Set a cutoff time: source the whole script when using a new cutoff value
cutoff <- 500
# Set the maximum time
max.time <- 365 * 3
# Parameter used for simulations 
target1.1.min <- 35
target1.1.max <- 50
target1.2.min <- 25
target1.2.max <- 50

target2.1.min <- 25
target2.1.max <- 30
target2.2.min <- 10
target2.2.max <- 30

target3.1.min <- 45
target3.1.max <- 80
target3.2.min <- 100
target3.2.max <- 100

niter <- 1000

K = 3
nc = c(25, 15, 10)
npt = 70
thetak = c(0.5, 5, 1)
rk = c(4, 1, 0.5, 3, 5, 1, 1, 1, 1)
lambda1 = 0.03
lambda2 = 0.04
max.time = 365 * 3
target.enroll.sub1 <- 250
target.enroll.sub2 <- 270


# generate data
# check there are enough enrolled patient until the end date
suppressWarnings(source('data_simu.R'))
n1 <- nrow(data1 %>% filter(status == 1 & screen.time <= max.time))
n2 <- nrow(data1 %>% filter(status == 2 & screen.time <= max.time))
while(n1 < target.enroll.sub1 | n2 < target.enroll.sub2){
  suppressWarnings(source('data_simu.R'))
  n1 <- nrow(data1 %>% filter(status == 1 & screen.time <= max.time))
  n2 <- nrow(data1 %>% filter(status == 2 & screen.time <= max.time))
}

## save data (before cutoff) for this seed
write.table(data1, paste0('../data/', 'data', seed, '.csv'), 
            col.names = TRUE, row.names = FALSE, sep = ",")


## Use the proposed Hierarchical Bayesian Method
## =====================================================================
## Step 2: source MCMC functions
## Below are main functions of our methods
source('Gibbs.R')
source('MH5.R')
source('label_update.R')
source('dynamic_fit_joint.R')

## =====================================================================
## Step 3: fit model
## set MCMC parameter
niter <- 3000
n.burnin <- 2001
set.seed(seed)
## run MCMC
## This steps takes a long time, around 10 mins.
iter.data <- jointfit(data1.final, niter)

## =====================================================================
## Step 4: prediction
## set prediction parameter
target.enroll <- target.enroll.sub1 + target.enroll.sub2
p.enroll.train <- 1- sum(data1.train$status == 3)/sum(data1.train$status %in% c(1, 2, 3))
target.screen <- target.enroll
target.screen.remain <- target.screen - sum(data1.train$screened == 1 & data1.train$status%in%c(1, 2))
n.repeat <- 1000
# only do prediction with centers:
# 1. planned before cutoff date
# 2. still open after cutoff date
# The above info is assumed to be known in reality
rule.1 <- center.info$center[center.info$plan == 1]
rule.2 <- center.info$center[center.info$close == 0]
center.include <- intersect(rule.1, rule.2)
J.pre <- length(center.include)
# sample h new IAT from each center
h <- round(target.screen.remain/J.pre * 150)
set.seed(seed)
suppressMessages(suppressWarnings(source('dynamic_predict_close.R')))

## =====================================================================
## Step 5: measurement
end.time <- max.time
source('measure.R')
# result with cap adjustment
write.table(matrix(c(seed, abs.dif1, res30$abs1.dif1, res60$abs1.dif1, res90$abs1.dif1), nrow = 1),
          '../result/result1.csv', col.names = FALSE, row.names = FALSE, sep = ',',
          append = TRUE)
# result without cap adjustment
write.table(matrix(c(seed, abs.dif, res30$abs.dif, res60$abs.dif, res90$abs.dif), nrow = 1),
          '../result/result0.csv', col.names = FALSE, row.names = FALSE, sep = ',',
          append = TRUE)

rm(iter.data, data.new, data.new.pre, data.new.pre0, data.new.pre.adj)

## After the completion of several repetitions, you can use result1.csv and 
## result0.csv to obtain results for BH and BHA in Table 2.
