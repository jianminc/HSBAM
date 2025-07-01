## required packages
library(dplyr)
library(coda)
library(MCMCpack)
library(survival)
library(mvtnorm)


## ggplot theme
myblue1 <- rgb(3/255, 154/255, 166/255)
myblue2 <- rgb(71/255, 92/255, 105/255)

mytheme <- theme_bw() +
  theme(
    # panel.border       = element_blank(),
    panel.grid.major.x = element_line(colour = "lightgray"),
    panel.grid.minor.x = element_line(colour = "lightgray"),
    panel.grid.major.y = element_line(colour = "lightgray"),
    panel.grid.minor.y = element_line(colour = "lightgray"),
    strip.background   = element_rect(fill = myblue2),
    strip.text         = element_text(colour = "white", face = "bold"),
    axis.title.y       = element_text(size = 10, margin = margin(0, 10, 0, 0)),
    axis.text.y        = element_text(size = 10),
    axis.text.x        = element_text(angle = 45, size = 10, hjust = 1),
    plot.title         = element_text(size = 11, hjust = 0.5, face = "bold"))

