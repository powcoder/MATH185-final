hubmean<- function(x,c)
{
  tol <- 10^(-6)
  n <- length(x)
  est_1 <- mean(x)
  est_2 <- 0
  iter <- 0
  
repeat{
     if(abs(est_2-est_1)<tol | iter>=200) break
      est_2 <- est_1
      r <- x - est_1
      a <- r*(abs(r)<=c) + sign(r)*c*(abs(r)>c)
      w <- a/r
      est_1 <- (w%*%x)/sum(w)
      iter <- iter + 1
}
  return(est_1)
}