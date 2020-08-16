library(truncnorm)
library(tidyverse)


# two sided selection, zcut is cutoff for |z|
OneSideConditionalTest = function(z, mu, zcut = 1.96, side = c("greater","less"), selectRegion = c("greater","both")) {
  side = match.arg(side)
  if(z<0) {
    z = -z
    mu = -mu
  }
  #if (z<zcut) zcut = z/2
  if (z<zcut) stop("|z| has to be more than the cutoff")
  selectregion = match.arg(selectRegion)
  if (selectregion == "both"){
    selectProb = 1 - pnorm(zcut, mu, 1) + pnorm(-zcut, mu, 1)
    if (side == "greater") {
      return(pnorm(z, mu, 1,lower.tail = FALSE)/selectProb)
    }
    else{
      return((pnorm(z,mu,1)-pnorm(zcut,mu,1)+pnorm(-zcut,mu,1))/selectProb)
    }
  } else 
  {
    if (side == "greater") {
      return(1 - ptruncnorm(z, zcut, Inf, mu, 1))
    }
    else{
      return(ptruncnorm(z,zcut, Inf, mu, 1))
    }
  }  
}


selectiveCI = function(delta, se, selectRegion = c("greater","both"), zcut = 1.96, cisize = 0.95) {
  alpha = 1 - cisize
  sn = sign(delta)
  z = abs(delta/se)
  if (z<zcut+0.2) zcut = z/2
  #if (z<zcut) stop("|z| has to be more than the cutoff")
  funGreater = function(x) OneSideConditionalTest(z, x, zcut,"greater",selectRegion) - alpha/2
  funLess = function(x) OneSideConditionalTest(z, x, zcut,"less",selectRegion) - alpha/2
  l = uniroot(funGreater, interval = c(-6.3, z))
  u = uniroot(funLess, interval = c(z, 1000))
  if(sn > 0){
    return(c(l$root, u$root)*se)
  } else {
    return(c(-u$root, -l$root)*se)
  }
}

#' delta must be positive and delta/se must be greater than or equal to zcut
selectiveMLE = function(delta, se, selectRegion = c("greater","both"), zcut = 1.96){
  sn = sign(delta)
  z = abs(delta/se)
  if (z<zcut) zcut = z/2
  #if (z<zcut) stop("|z| has to be more than the cutoff")
  selectregion = match.arg(selectRegion)
  if (selectregion == "both")
  {
    func = function(mu){
      return(z - mu - (dnorm(zcut-mu)-dnorm(-zcut-mu))/(pnorm(zcut-mu,lower.tail = FALSE)+pnorm(-zcut-mu)))
    }
    mle = uniroot(func, interval = c(-10,2*z))
  } else {
    func = function(mu){
      return(z - mu - (dnorm(zcut-mu))/(pnorm(zcut-mu,lower.tail = FALSE)))
    }
    mle = uniroot(func, interval = c(-30,2*z))
  }
  return(sn* mle$root*se)
}
