hubreg<-function(X, Y, c)
# Use Huber regression to fit a linear model Y = beta0 + X*beta + epsilon 
# X: n by d design matrix
# Y: n by 1 response vector
# beta0: intercept
# beta: d by 1 vector of regression coefficients  
# c: tuning parameter in the Huber loss  

## Output: 
# d+1 by 1 vector of estimated beta0 and beta
# number of iterations  
    
{
  tol<-10^(-6)
  n<-nrow(X)  
  d<-ncol(X)
  beta_old<-matrix(0, nr=d+1, nc=1)
  beta_new<-lm(Y~X)$coef  
  Z<-cbind(matrix(1, nr=n, nc=1), X)
  res<-Y-Z%*%beta_new
  
  w_old<-matrix(0, nr=n, nc=1)
  w_new<-matrix(1, nr=n, nc=1)
  iter<-0
  
  repeat{
    if( max(abs(beta_old - beta_new)) < tol | iter>200) break
    beta_old<-beta_new
    w_old<-w_new
    a<-rep(0, n)
    a[which(abs(res)>=c)]<-1
    w_new<-(abs(res)/c - 1)*a
    D<-diag(c(1/(1+w_new)))
    beta_new<-solve(t(Z)%*%D%*%Z/n)%*%(t(Z)%*%D%*%Y/n)
    res<-Y-Z%*%beta_new
    iter<-iter+1
  }
  
  beta<-matrix(beta_new)
  outlist<-list(beta=beta, iter=iter)
  return(outlist)
}