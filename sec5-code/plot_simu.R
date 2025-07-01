# Visualization for the true data
#========================================================================
## Overall cumulative number of screened patients vs time
## Training
plotmat <- data1.final %>% filter(screened == 1) %>% 
  dplyr::select(screened, screen.time) %>% 
  arrange(screen.time)
plotmat$cum.screened <- cumsum(plotmat$screened)

p.screen.train <- ggplot(plotmat, aes(screen.time, cum.screened)) + 
  geom_line() + 
  xlim(c(0, max(plotmat$screen.time) + 5)) + 
  ylab(c('Number of screened patients')) + 
  xlab(c('Time')) + 
  ggtitle('Training set: Overall Screening') + 
  mytheme


## Add testing data
plotmat$source <- 'train'
plotmat2 <- data1.test %>% 
  dplyr::select(screened, screen.time) %>% arrange(screen.time)
plotmat2$screened <- 1
plotmat2$cum.screened <- cumsum(plotmat2$screened) + max(plotmat$cum.screened)
plotmat2$source <- 'test'
plotmat2 <- rbind(plotmat, plotmat2)
p.screen.all <- ggplot(plotmat2, aes(screen.time, cum.screened, col = source)) + 
  geom_line()+ 
  xlim(c(0, max(plotmat2$screen.time) + 5))+ 
  ylab(c('Number of screened patients')) + 
  xlab(c('Time')) + 
  ggtitle('Whole set: Overall Screening') + 
  mytheme

##========================================================================
plotmat.screen.enroll <- data1.final %>% filter(screened == 1) %>% 
  filter(status %in% c(1, 2)) %>% 
  dplyr::select(screened, screen.time) %>% 
  arrange(screen.time)
plotmat.screen.enroll$cum.enroll <- cumsum(plotmat.screen.enroll$screened)
p.screen.enroll.train <- ggplot(plotmat.screen.enroll, aes(screen.time, cum.enroll)) + 
  geom_line() + 
  xlim(c(0, max(plotmat.screen.enroll$screen.time) + 5)) +
  ylab(c('Number of enrolled patients')) + 
  xlab(c('Time')) + 
  ggtitle('Training set: Overall Enrollment') + 
  mytheme


## add data for test
plotmat.screen.enroll$source <- 'train'
plotmat.screen.enroll2 <- data1.test %>% 
  filter(status %in% c(1, 2)) %>% 
  dplyr::select(screened, screen.time) %>% 
  arrange(screen.time)
plotmat.screen.enroll2$screened <- 1
plotmat.screen.enroll2$cum.enroll <- cumsum(plotmat.screen.enroll2$screened) + max(plotmat.screen.enroll$cum.enroll)
plotmat.screen.enroll2$source <- 'test'
plotmat.screen.enroll2 <- rbind(plotmat.screen.enroll, plotmat.screen.enroll2)
p.screen.enroll.all <- ggplot(plotmat.screen.enroll2, 
                              aes(screen.time, cum.enroll, col = source)) + 
  geom_line() + 
  xlim(c(0, max(plotmat.screen.enroll2$screen.time) + 5)) +
  ylab(c('Number of enrolled patients')) + 
  xlab(c('Time')) + 
  ggtitle('Whole set: Overall Enrollment') + 
  mytheme
##========================================================================
## Overall cumulative number of enrolled patients vs time
plotmat.enroll <- data1.final %>% filter(screened == 1) %>% 
  filter(status %in% c(1, 2)) %>% 
  dplyr::select(screened, enroll.time) %>% 
  arrange(enroll.time)
plotmat.enroll$cum.enroll <- cumsum(plotmat.enroll$screened)
p.enroll.train <- ggplot(plotmat.enroll, aes(enroll.time, cum.enroll)) + 
  geom_line() + 
  xlim(c(0, max(plotmat.enroll$enroll.time) + 5)) +
  ylab(c('Number of enrolled patients')) + 
  xlab(c('Time')) + 
  ggtitle('Training set: Overall Enrollment') + 
  mytheme


## add data for test
plotmat.enroll$source <- 'train'
plotmat.enroll2 <- data1.test %>% 
  filter(status %in% c(1, 2)) %>% 
  dplyr::select(screened, enroll.time) %>% 
  arrange(enroll.time)
plotmat.enroll2$screened <- 1
plotmat.enroll2$cum.enroll <- cumsum(plotmat.enroll2$screened) + max(plotmat.enroll$cum.enroll)
plotmat.enroll2$source <- 'test'
plotmat.enroll2 <- rbind(plotmat.enroll, plotmat.enroll2)
p.enroll.all <- ggplot(plotmat.enroll2, 
                       aes(enroll.time, cum.enroll, col = source)) + 
  geom_line() + 
  xlim(c(0, max(plotmat.enroll2$enroll.time) + 5)) +
  ylab(c('Number of enrolled patients')) + 
  xlab(c('Time')) + 
  ggtitle('Whole set: Overall Enrollment') + 
  mytheme


## ========================================================================
## look at enrollment for each subtype
## type 1
plotmat.enroll.type1 <- data1.final %>% filter(screened == 1) %>% 
  filter(status == 1) %>% 
  dplyr::select(screened, enroll.time) %>% 
  arrange(enroll.time)
plotmat.enroll.type1$cum.enroll <- cumsum(plotmat.enroll.type1$screened)
p.enroll.train.st1 <- ggplot(plotmat.enroll.type1, 
                             aes(enroll.time, cum.enroll)) + 
  geom_line() + 
  xlim(c(0, max(plotmat.enroll.type1$enroll.time) + 5)) +
  ylab(c('Number of enrolled patients')) + 
  xlab(c('Time')) + 
  ggtitle('Training set: Sub-type1 Enrollment') + 
  mytheme

## add data for test
plotmat.enroll.type1$source <- 'train'
plotmat.enroll2.type1 <- data1.test %>% 
  filter(status == 1) %>% 
  dplyr::select(screened, enroll.time) %>% 
  arrange(enroll.time)
plotmat.enroll2.type1$screened <- 1
plotmat.enroll2.type1$cum.enroll <- cumsum(plotmat.enroll2.type1$screened) + max(plotmat.enroll.type1$cum.enroll)
plotmat.enroll2.type1$source <- 'test'
plotmat.enroll2.type1 <- rbind(plotmat.enroll.type1, plotmat.enroll2.type1)
p.enroll.all.st1 <- ggplot(plotmat.enroll2.type1, 
                           aes(enroll.time, cum.enroll, col = source)) + 
  geom_line() + 
  xlim(c(0, max(plotmat.enroll2.type1$enroll.time) + 5)) +
  ylab(c('Number of enrolled patients')) + 
  xlab(c('Time')) + 
  ggtitle('Whole set: Sub-type1 Enrollment') + 
  mytheme

##==========================================================================
## type 2
plotmat.enroll.type2 <- data1.final %>% filter(screened == 1) %>% 
  filter(status == 2) %>% 
  dplyr::select(screened, enroll.time) %>% 
  arrange(enroll.time)
plotmat.enroll.type2$cum.enroll <- cumsum(plotmat.enroll.type2$screened)
p.enroll.train.st2 <- ggplot(plotmat.enroll.type2, aes(enroll.time, cum.enroll)) + 
  geom_line() + 
  xlim(c(0, max(plotmat.enroll.type2$enroll.time) + 5)) +
  ylab(c('Number of enrolled patients')) + 
  xlab(c('Time')) + 
  ggtitle('Training set: Sub-type2 Enrollment') + 
  mytheme

## add data for test
plotmat.enroll.type2$source <- 'train'
plotmat.enroll2.type2 <- data1.test %>% 
  filter(status == 2) %>% 
  dplyr::select(screened, enroll.time) %>% 
  arrange(enroll.time)
plotmat.enroll2.type2$screened <- 1
plotmat.enroll2.type2$cum.enroll <- cumsum(plotmat.enroll2.type2$screened) + max(plotmat.enroll.type2$cum.enroll)
plotmat.enroll2.type2$source <- 'test'
plotmat.enroll2.type2 <- rbind(plotmat.enroll.type2, plotmat.enroll2.type2)
p.enroll.all.st2 <- ggplot(plotmat.enroll2.type2, aes(enroll.time, cum.enroll, col = source)) + 
  geom_line() + 
  xlim(c(0, max(plotmat.enroll2.type2$enroll.time) + 5)) +
  ylab(c('Number of enrolled patients')) + 
  xlab(c('Time')) + 
  ggtitle('Whole set: Sub-type2 Enrollment') + 
  mytheme

