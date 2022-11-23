lpreg<-function(x,y,x0,deg,h)
# Use local polynomial regression to fit a model y = m(x) + epsilon 
# (x,y): input observations
# x0: the point at which we want to estimate -- m(x0)
# deg: degree of the polynomial (deg=1:linear;deg=2:quadratic;deg=3:cubic)
# h: bandwidth
# Output: estimated value of m(x0)
{ 
  n <- length(x)  # sample size
  if(x0<min(x)-h | x0>max(x)+h){
     lpfit <- 0
  } else{
     X <- matrix(1, nr=n, nc=1)
    for(j in 1:deg){
     X <- cbind(X, (x-x0)^j/factorial(j))
    }                              # n by (deg+1) design matrix
     W <- diag(dnorm((x-x0)/h))    # n by n diagonal matrix
     coef <- solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%y # solve weighted least squares
     lpfit <- coef[1]
  }
  return(lpfit)
}