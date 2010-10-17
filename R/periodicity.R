#
# periodicity.R
#
# copyright (c) 2010, Danny Arends
# last modified Okt, 2010
# first written Okt, 2010
# 
# R functions: findperiod
# Basic functions to explore periodic phenomena
#

findperiod <- function(x, y, N = max(x)-min(x)){
  fftx <- fft(y)
  periodogram <- Mod(fftx)^2/N
  periodogram <- periodogram[-1]
  periodogram <- periodogram[1:(length(y)/2)]
  op <- par(mfrow = c(2,1))
  plot(x,y,type="l")
  plot(N/1:(length(y)/2),periodogram)
  r <-  N / which.max(periodogram)
  r
}

#findperiod(seq(0,10,0.1),sin(3+seq(0,10,0.1)))