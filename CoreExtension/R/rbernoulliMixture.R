rbernoulliMixture <-
function(pi,theta,n){
  ## Simulations of observations x
  ## according to mixture of r different size-p multivariate Bernoulli

  r<-length(pi)  # number of component	
  p<-ncol(theta) # dimension of the observations

  ## Simulate the component origins
  z<-rmultinom(1, size=n, prob = pi)
  z<-rep(1:r,z)
  
  ##Simulate the data matrix
  x<-matrix(0,n,p) # data matrix
  for (i in 1:n)
    for(j in 1:p)
      x[i,j]<-rbinom(1,1,theta[z[i],j])

  return(list(x=x,z=z,pi=pi,theta=theta))
}
