remove(list = ls())
library(KernSmooth)
library(NonpModelCheck)
source("lpreg.r")
cmb <- read.csv('cmb.csv')
freq <- cmb$freq
power <- cmb$power
plot(freq,power)
h1 <- dpill(freq,power)


grid <- seq(min(freq),max(freq),1.5)
M <- length(grid)

## Local liear regression
fit1 <- numeric(M)
for (m in 1:M){
  fit1[m] <- lpreg(freq,power,grid[m],1,h1)
}
lines(grid,fit1,col='green',lwd=2)

## Local cubic regression
fit3 <- numeric(M)
for (m in 1:M){
  fit3[m] <- lpreg(freq,power,grid[m],3,h1)
}
lines(grid,fit3,col='red',lwd=2)


## Preprocessing
x <- unique(freq)
n <- length(x)
y <- numeric(n)
for(i in 1:n){
   y[i] <- mean(power[freq==x[i]])
}
plot(x,y)

## Compute LOO-MSE for NW estimator
band_grid <- seq(1,100,0.5)
L <- length(band_grid)
get_nwmse <- function(x,y,h)
{
  mse <- 0
  n <- length(x)
  for(i in 1:n){
    xi <- x[-i]   # remove i-th obs
    yi <- y[-i]   # remove i-th obs
    fiti <- ksmooth(xi,yi,kernel='normal',bandwidth=h,
                    x.points=x[i]) # compute the fit at i-th obs
    mse <- mse + (fiti$y - y[i])^2
  }
  return(mse/n)
}
nwmse <- numeric(L)
for (k in 1:L){
  nwmse[k] = get_nwmse(x,y,band_grid[k])
}
plot(nwmse,type='l',col='red', lwd=2)
nw_h <- band_grid[which(nwmse == min(nwmse))]
nw_fit <- ksmooth(x,y,kernel='normal',bandwidth=2*nw_h)
plot(x,y)
lines(nw_fit,col='blue',lwd=2)


## Compute LOO-MSE for local linear estimator
band_grid <- seq(1,100,0.5)
L <- length(band_grid)
get_llmse <- function(x,y,h)
{
  mse <- 0
  n <- length(x)
  for(i in 1:n){
    xi <- x[-i]
    yi <- y[-i]  
    fiti <- lpreg(xi,yi,x[i],1,h)
    mse <- mse + (fiti-y[i])^2
  }
  return(mse/n)
}
llmse <- numeric(L)
for (k in 1:L){
  llmse[k] = get_llmse(x,y,band_grid[k])
}
plot(llmse,type='l',col='red', lwd=2)
ll_h <- band_grid[which(llmse == min(llmse))]

h1 <- dpill(x,y)
ll_fit <- numeric(n)
for (i in 1:n){
  ll_fit[i] <- lpreg(x,y,x[i],1,ll_h)
}
plot(x,y)
lines(x,ll_fit,col='green',lwd=2)
nw_fit <- ksmooth(x,y,kernel='normal',bandwidth=2*h1)
lines(nw_fit,col='red',lwd=2)

