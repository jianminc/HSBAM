rm(list = setdiff(ls(), c()))
source('preamble.R')
## =============================================================================
## Step 1: generate simulated data
set.seed(5)
# Set a cutoff time: source the whole script when using a new cutoff value
cutoff <- 300
# Set the maximum time
max.time <- 365 * 3
# 10 centers in region1
target1.1.min <- 35
target1.1.max <- 50
target1.2.min <- 25
target1.2.max <- 50
# 5 centers in region2
target2.1.min <- 25
target2.1.max <- 30
target2.2.min <- 15
target2.2.max <- 30
# 15 centers in region3
target3.1.min <- 70
target3.1.max <- 100
target3.2.min <- 100
target3.2.max <- 120
# 20 centers in region4
target4.1.max <- 150
target4.2.max <- 150
# Parameter used for simulations are be store in the following sourced file
source('simucode.R')
J # exact number of centers in the simulated data
# Genevrate plot to see the format of the simulated data
source('plot_simu.R')
# overall screening
p.screen.all
p.enroll.all
p.enroll.all.st1
p.enroll.all.st2
p.screen.enroll.all
grid.arrange(p.screen.train, p.screen.all, nrow = 2)
# overall enrollment
grid.arrange(p.enroll.train, p.enroll.all, nrow = 2)
# subtype enrollment
grid.arrange(p.enroll.all.st1, p.enroll.all.st2, nrow = 2)

## =====================================================================
## Step 2: source MCMC functions
source('Gibbs.R')
source('AM.R')
source('label_update.R')

source('Gibbs.R')
source('MH5.R')
source('label_update.R')
source('dynamic_fit_joint.R')

## =====================================================================
## Step 3: fit model
## set MCMC parameter
niter <- 5000
n.burnin <- 2001
set.seed(1)
source('dynamic_fit_joint.R')
iter.data <- jointfit(data1.final, niter = niter)

## =====================================================================
## Step 4: check convergence
plot(n.burnin:niter, log(iter.data$theta[n.burnin:niter, 2]), type='l')
plot(n.burnin:niter, log(iter.data$gamma[1, 1, n.burnin:niter]), type='l')
plot(n.burnin:niter, iter.data$p[10, 1, n.burnin:niter], type='l')
plot(n.burnin:niter, iter.data$u[n.burnin:niter, 1], type='l')
plot(n.burnin:niter, iter.data$lambda[n.burnin:niter], type='l')
## =====================================================================
## Step 5: prediction
## set prediction parameter
target.enroll.sub1 <- 230
target.enroll.sub2 <- 230
target.enroll <- target.enroll.sub1 + target.enroll.sub2
p.enroll.train <- 1- sum(data1.train$status == 3)/sum(data1.train$status %in% c(1, 2, 3))
target.screen <- target.enroll / p.enroll.train
target.screen.remain <- target.screen - sum(data1.train$screened == 1)
n.repeat <- 100
# only do prediction with centers: 
# 1. planned before cutoff date
# 2. still open after cutoff date
# The above info is assumed to be known in reality
rule.1 <- center.info$center[center.info$plan == 1]
rule.2 <- center.info$center[center.info$close == 0]
center.include <- intersect(rule.1, rule.2)
J.pre <- length(center.include)
# sample h new IAT from each center
h <- round(target.screen.remain/J.pre * 10)
set.seed(1)
#source('dynamic_predict_new.R')
source('dynamic_predict_close.R')

## =====================================================================
## Step 6: result visualization
end.time <- 900
source('result.R')
pre.screen
prob.screen
pre.enroll
prob.enroll
pre.enroll.st1
prob.enroll.st1
pre.enroll.st2
prob.enroll.st2

