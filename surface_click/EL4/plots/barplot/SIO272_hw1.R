library(magrittr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(reshape2)
library(pwr)

setwd("/home/sherry/SIO272")
data <- read.csv("bohar_2014-SIO272_2019.csv")

#--------Q1----------
meanSize_2014 <- mean(data$length_bohar_2014_mm)
cat("The mean size of Lutjanus bohar in 2014 is",meanSize_2014,"mm ,which is",meanSize_2014-380,"mm bigger than size of maturity (380mm)")

#Question2
library("pwr")
sampleSize <- nrow(data)
#power = 1-beta(type2 error)
d_val <- 20/ars(data$sampleSize)
pwr.t.test(d = 20,sig.level = 0.05,power = 0.9,type ="one.sample")

