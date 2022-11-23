# Math 185 HW1 WIN2019

### Problem 5 ###

unexposed=c(8,11,12,14,20,43,111)
exposed=c(35,56,83,92,128,150,176,208)

### option I: self-coding:

# tests if H0: x is distributed same as y
# H1:x is stochastically smaller then y
wilcoxwrite = function(x,y){
  count=0
  n=length(x)
  m=length(y)
  for (i in 1:n){
    for (j in 1:m){
      if (x[i]<y[j]){count=count+1}      
    }
  }
  # this uses the equality of Wilcoxin and Mann-Whitney U
  result=count+0.5*m*(m+1)
  return(result)
}

wilcoxwrite((unexposed+25), exposed)
# [1] 81
# notice canceling the 0.5*m*(m+1) will also give 45(as in later built-in solution)

# apply permutation to find approximation of p-value
# H0: exposed mean is same as unexposed+25 mean
# H1: unexposed+25 mean is less then exposed mean
B=100000
pvalwrite=function(xx,yy){
  n=length(xx)
  m=length(yy)
  count=0
  w0=wilcoxwrite(xx,yy)

  for (b in 1:B){
    z=sample(c(xx,yy), n+m , replace=FALSE)
    newx=z[1:n]
    newy=z[(n+1):  (n+m)]
    w=wilcoxwrite(newx, newy)
    if (w>=w0){count=count+1}    
  }
  return((count+1)/(B+1))
}
une25=unexposed+25
pvalwrite(une25, exposed)
# [1]  0.02745973
# reject H0, and conclude at 5% that exposed mean is higher then unexposed+25 mean

### option II: built in function:
# Note the below is actually the Mann Whitney U Test
wilcox.test(exposed, (unexposed+25) , alternative = "g") 
# It will print the following:
# Wilcoxon rank sum test
# data:  exposed and (unexposed + 25)
# W = 45, p-value = 0.02704
# alternative hypothesis: true location shift is greater than 0
# also reject H0 at 5% level

### to form 95% CI for the difference delta, we can use bootstrap
B=100000
diff=numeric(B)
for (b in 1:B){
  x=sample(exposed, length(exposed), replace=TRUE)
  y=sample(unexposed, length(unexposed), replace=TRUE)
  diff[b]=mean(x)-mean(y)  
}
quantile(diff, c(0.025, 0.0975)) 
# thus the 95% confidence interval for difference in mean is
# 2.5%    9.75% 
# 37.53571 53.73214 
# notice 25 should not be inside


