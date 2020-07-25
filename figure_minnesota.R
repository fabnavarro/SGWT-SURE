rm(list=ls())
library(gasper)
library(ggplot2)
library(gridExtra)

#- Load signals
load("minnesota/fminnesota.Rdata")
f1 <- fminnesota$f1
f2 <- fminnesota$f2

#- Graph of the minnesota with the signals added
p1 <- plot_signal(minnesota, f1, size=15*f1)
p2 <- plot_signal(minnesota, f2, size=15*f2)
#grid.arrange(p1, p2, ncol=2)
