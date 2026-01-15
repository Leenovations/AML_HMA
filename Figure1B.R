library(ggplot2)
library(tidyverse)
library(viridis)
library(plotly)
library(manipulate)
library(openair)
library(RColorBrewer)
#------------------------------------------------------------------------------#
Data <- read.table('Figure1B.txt',header=T)
Data <- na.omit(Data)
#------------------------------------------------------------------------------#
viridis <- c('#440154', '#3b528b', '#21918c', '#5ec962', '#fde725')
#------------------------------------------------------------------------------#
smoothScatter(Data$CR, Data$NR,
              transformation = function(x) x ^0.25,
              colramp = colorRampPalette(c(viridis)),
              cex= 0.75,
              pch=NA,
              xlab='CR',
              ylab='NR',
              bandwidth = 0.05)
#------------------------------------------------------------------------------#