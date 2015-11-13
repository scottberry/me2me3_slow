rm(list=ls())

## test data for different noise inputs
noise = c(1,10,100,200,500,1000)
SD = c(0.04,0.09,0.28,0.4,0.6,0.9)

## Plot log-log
par(mfrow=c(1,1),mar=c(2,4,0,0)+0.5,oma=c(3,3,0,0))
plot(noise,SD,log="xy")

lm(log(SD) ~ log(noise))

# log(SD) = 0.4539*log(noise) - 3.3215
# SD = noise^0.4539 + 0.03609864
SDpred <- function(x) {
  exp(-3.3215)*x^0.4539 
}

x <- seq(1,1000)
y <- sapply(x,FUN = SDpred)

plot(noise,SD)
lines(x,y)

## Invert function
NoiseReq <- function(y) {
  (y/exp(-3.3215))^(1/0.4539)
}

y <- seq(0,1,by=0.05)
x <- sapply(y,FUN = NoiseReq)

plot(SD,noise)
lines(y,x)

y <- seq(0,1.0,by=0.05)
x <- sapply(y,FUN = NoiseReq)
x

