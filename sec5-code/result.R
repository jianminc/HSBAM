## make plot for prediction
data.new.pre <- data.new.pre.adj

# ============================================================
# enrollment for sub-type 1
date_min <- min(sapply(data.new.pre, FUN = function(x) {
  min(x$enroll.time[x$status == 1], na.rm = TRUE)}))
date_max <- max(sapply(data.new.pre, FUN = function(x) {
  max(x$enroll.time[x$status == 1], na.rm = TRUE)}))
dates <- seq(date_min, date_max, by = 1/30)
dates <- unique(c(dates, date_max))
# NUMBER OF ARRIVALS WITHIN EACH INTERVAL
counts <- sapply(data.new.pre, FUN = function(x)
  hist(x      = x$enroll.time[x$status == 1],
       plot   = FALSE,
       breaks = dates,
       right  = FALSE)$counts
)

# CUMULATIVE NUMBER OF ARRIVALS OVER TIME
counts_cum <- apply(counts, MARGIN = 2, FUN = cumsum)
# SUMMARISE (MEDIAN AND QUANTILES)
df <- data.frame(date   = dates[-1],
                 median = apply(counts_cum, MARGIN = 1, FUN = median),
                 upp    = apply(counts_cum, MARGIN = 1, FUN = quantile, 0.975),
                 low    = apply(counts_cum, MARGIN = 1, FUN = quantile, 0.025),
                 upp80    = apply(counts_cum, MARGIN = 1, FUN = quantile, 0.10),
                 low80    = apply(counts_cum, MARGIN = 1, FUN = quantile, 0.90),
                 upp60    = apply(counts_cum, MARGIN = 1, FUN = quantile, 0.20),
                 low60    = apply(counts_cum, MARGIN = 1, FUN = quantile, 0.80))
df$date <- df$date/30

## get the 95% predictive interval 
complete.time.pre <- rep(0, n.repeat)
for(i in 1:n.repeat){
  pre <- data.new.pre[[i]] %>% filter(status == 1) %>% 
    arrange(enroll.time) %>% slice(target.enroll.sub1:target.enroll.sub1) %>% 
    dplyr::select(enroll.time)
  complete.time.pre[i] <- unlist(pre)
}
q975 <- round(quantile(complete.time.pre, 0.975), 2) /30
q025 <- round(quantile(complete.time.pre, 0.025), 2) /30

cutoffn <- max(plotmat.enroll2.type1 %>% filter(enroll.time <= cutoff) %>%
  dplyr::select(cum.enroll))
ggplot(data = df, aes(x = date, y = median)) +
  geom_ribbon(aes(ymin = low, ymax = upp, fill = "95% prediction interval"), alpha = 0.6) +
  geom_ribbon(aes(ymin = low80, ymax = upp80, fill = "80% prediction interval"), alpha = 0) +
  geom_ribbon(aes(ymin = low60, ymax = upp60, fill = "60% prediction interval"), alpha = 0) +
  geom_line(colour = 'blue', size = 0.7) +
  geom_vline(xintercept = as.numeric(cutoff)/30, colour = "darkgray",
             size = 0.7, linetype = "longdash") + 
  geom_line(data = plotmat.enroll2.type1[1:target.enroll.sub1, ], aes(x = enroll.time/30, y = cum.enroll), 
            color = 'red'
            , size = 0.7) + 
  xlim(c(0, 35)) + 
  ylim(c(0, 450)) +
  ylab(c('Number of enrolled patients')) + 
  xlab(c('Time (months)')) + 
  ggtitle(paste0(round(cutoff/30),' months: Subtype1 Enrollment')) + 
  mytheme + 
  theme(legend.position="none") + 
  annotate(geom = "text", x = plotmat.enroll2.type1$enroll.time[plotmat.enroll2.type1$cum.enroll == target.enroll.sub1]/30,
           y = target.enroll.sub1, 
           label = "x", size = 5, colour = "black") +
  annotate("text",
           x      = cutoff/30,
           y      = 0,
           hjust  = 0,
           label  = paste0('cutoff: ', cutoffn, 'pts/', floor(cutoffn/target.enroll.sub1*100), '%'),
           parse  = FALSE,
           size   = 3,
           colour = "black")+
  annotate("text",
           x      = 0,
           y      = target.enroll.sub1 + 30,
           hjust  = 0,
           label  = paste0('Target:', target.enroll.sub1),
           parse  = TRUE,
           size   = 3,
           colour = "red")+
  geom_hline(yintercept = target.enroll.sub1, color = 'red', 
             size = 0.3, linetype = "longdash")



# ============================================================
# enrollment for sub-type 2
date_min <- min(sapply(data.new.pre, FUN = function(x) {
  min(x$enroll.time[x$status == 2], na.rm = TRUE)}))
date_max <- max(sapply(data.new.pre, FUN = function(x) {
  max(x$enroll.time[x$status == 2], na.rm = TRUE)}))
dates <- seq(date_min, date_max, by = 1/30)
dates <- unique(c(dates, date_max))
# NUMBER OF ARRIVALS WITHIN EACH INTERVAL
counts <- sapply(data.new.pre, FUN = function(x)
  hist(x      = x$enroll.time[x$status == 2],
       plot   = FALSE,
       breaks = dates,
       right  = FALSE)$counts
)

# CUMULATIVE NUMBER OF ARRIVALS OVER TIME
counts_cum <- apply(counts, MARGIN = 2, FUN = cumsum)
# SUMMARISE (MEDIAN AND QUANTILES)
df <- data.frame(date   = dates[-1],
                 median = apply(counts_cum, MARGIN = 1, FUN = median),
                 upp    = apply(counts_cum, MARGIN = 1, FUN = quantile, 0.975),
                 low    = apply(counts_cum, MARGIN = 1, FUN = quantile, 0.025),
                 upp80    = apply(counts_cum, MARGIN = 1, FUN = quantile, 0.10),
                 low80    = apply(counts_cum, MARGIN = 1, FUN = quantile, 0.90),
                 upp60    = apply(counts_cum, MARGIN = 1, FUN = quantile, 0.20),
                 low60    = apply(counts_cum, MARGIN = 1, FUN = quantile, 0.80))
df$date <- df$date/30

## get the 95% predictive interval 
complete.time.pre <- rep(0, n.repeat)
for(i in 1:n.repeat){
  pre <- data.new.pre[[i]] %>% filter(status == 2) %>% 
    arrange(enroll.time) %>% slice(target.enroll.sub2:target.enroll.sub2) %>% 
    dplyr::select(enroll.time)
  complete.time.pre[i] <- unlist(pre)
}
q975 <- round(quantile(complete.time.pre, 0.975), 2) /30
q025 <- round(quantile(complete.time.pre, 0.025), 2) /30

cutoffn <- max(plotmat.enroll2.type2 %>% filter(enroll.time <= cutoff) %>%
                 dplyr::select(cum.enroll))
ggplot(data = df, aes(x = date, y = median)) +
  geom_ribbon(aes(ymin = low, ymax = upp, fill = "95% prediction interval"), alpha = 0.6) +
  geom_ribbon(aes(ymin = low80, ymax = upp80, fill = "80% prediction interval"), alpha = 0) +
  geom_ribbon(aes(ymin = low60, ymax = upp60, fill = "60% prediction interval"), alpha = 0) +
  geom_line(colour = 'blue', size = 0.7) +
  geom_vline(xintercept = as.numeric(cutoff)/30, colour = "darkgray",
             size = 0.7, linetype = "longdash") + 
  geom_line(data = plotmat.enroll2.type2[1:target.enroll.sub2, ], aes(x = enroll.time/30, y = cum.enroll), 
            color = 'red'
            , size = 0.7) + 
  xlim(c(0, 35)) + 
  ylim(c(0, 450)) +
  ylab(c('Number of enrolled patients')) + 
  xlab(c('Time (months)')) + 
  ggtitle(paste0(round(cutoff/30),' months: Subtype2 Enrollment')) + 
  mytheme + 
  theme(legend.position="none") + 
  annotate(geom = "text", 
           x = plotmat.enroll2.type2$enroll.time[plotmat.enroll2.type2$cum.enroll == target.enroll.sub2/30],
           y = target.enroll.sub2, 
           label = "x", size = 5, colour = "black") +
  annotate("text",
           x      = cutoff/30,
           y      = 0,
           hjust  = 0,
           label  = paste0('cutoff: ', cutoffn, 'pts/', floor(cutoffn/target.enroll.sub2*100), '%'),
           parse  = FALSE,
           size   = 3,
           colour = "black") +
  annotate("text",
           x      = 0,
           y      = target.enroll.sub2 + 30,
           hjust  = 0,
           label  = paste0('Target:', target.enroll.sub2),
           parse  = TRUE,
           size   = 3,
           colour = "red")+
  geom_hline(yintercept = target.enroll.sub2, color = 'red', 
             size = 0.3, linetype = "longdash")



