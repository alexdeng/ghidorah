#required packages: tidyverse, R6, caret, xgboost, truncnorm, rmutil, Rfast

library(tidyverse)


##Simulation Study##

#'add some proportion of 0s
r_zeroInflated = function(n, p, generator){
  nzero = floor(n*p)
  x = generator(n-nzero)
  #c(rep(0,nzero),x)
  sample(c(rep(0,nzero),x))
}

#' simulate random t-distribution, sd here is the standard deviation (different from scale parameter of t)
rt_scale = function(n, df, ncp, sd){
  sd/sqrt(df/(df-2))*rt(n, df, ncp)
}

#' simulate random laplace distribution, sd here is the standard deviation
rlaplace_scale = function(n, sd){
  rmutil::rlaplace(n, 0, sd/sqrt(2))
}

#'simulate random Huber distribution, sd 
rhuber <- function(n, mu, b, K){
  sigma = sqrt(K/b)
  
  C = sqrt(2*pi*sigma^2)*(pnorm(K/sigma)- pnorm(-K/sigma)) + 2*sigma^2/K*exp(-K^2/2/sigma^2)
  C1 = exp(-K^2/2/sigma^2)*sigma^2/K/C
  C2 = C1 + sqrt(2*pi*sigma^2)*(pnorm(K/sigma)- pnorm(-K/sigma))/C
  U = runif(n)
  ind1= U <C1
  ind2 = C1<U & U<C2
  ind3 = U > C2
  y = rep(0, n)
  y[ind1] = (sigma^2/K*log(U[ind1]*C*K/sigma^2)-K/2+mu)
  y[ind2] = qnorm((C*U[ind2]-exp(-K^2/2/sigma^2)*sigma^2/K)/sqrt(2*pi*sigma^2) + pnorm(-K/sigma))*sigma + mu
  y[ind3] = (log( (-C*U[ind3] + sqrt(2*pi*sigma^2)*(pnorm(K/sigma)- pnorm(-K/sigma)) +  exp(-K^2/2/sigma^2)*sigma^2/K)*K/sigma^2 + exp(-K^2/2/sigma^2) )-K^2/2/sigma^2)*(-sigma^2/K)+mu  
  return(y)
}


#' simulate data used for comparison
#' @param m Number of data points (experiment) to simulate
#' @param priorGen a function takes m as parametr and simulate m effect sizes
#' @param sampleSizeGen a function takes m as pamater and simulate m sample sizes. Scale^2/sampleSize is the variance of the noise. 
#' @param scale the standard deviation of the noise before divided by sample size
#' @return a tibble with the following schema: y(observation), mu(true effect),y_A, y_B (splitted observations), y_scaled, y_scaled_A, y_scaled_B, mu_scaled (observations scaled by pooled standard deviation, mu_scaled is the effect size), sample_size (efficient sample size N), sample_size_A, sample_size_B (sample size of two splits), scale (the pooled standard deviation), scale_A, scale_B. 
#' Scale_A and scale_B are the same as scale in this case. In empirical data they might be different. 
#' The prior is defined on top of effect size --- effect divided by the pooled standard deviation (scale). When divided by the scale, noises have variance 1/sample_size.
simulateData = function(m, priorGen, sampleSizeGen, scale){
  mu_scaled = priorGen(m)
  n = sampleSizeGen(m)
  # noises for two equal splits
  eps_A = rnorm(m, 0, sqrt(2/n))
  eps_B = rnorm(m, 0, sqrt(2/n))
  eps = (eps_A + eps_B)/2
  y_scaled_A = mu_scaled + eps_A
  y_scaled_B = mu_scaled + eps_B
  y_scaled = mu_scaled + eps
  y_A = y_scaled_A * scale
  y_B = y_scaled_B * scale
  y = y_scaled * scale
  mu = mu_scaled * scale
  # tstat = y_scaled*sqrt(n)
  # pval = 2*pnorm(abs(tstat),lower=FALSE)
  tibble::tibble(mu, y, y_A, y_B, y_scaled, y_scaled_A, y_scaled_B, mu_scaled, sample_size = n, sample_size_A = n/2, sample_size_B = n/2, scale, scale_A=scale, scale_B = scale)
}


#' Main evaluation entry function to evaluate based on RMSE
#' @param methods is a list of R6 classes for different techniques/methods to be evaluated
evaluateSim = function(B, trainSize, testSize, priorGen, sampleSizeGen, scale, methods, trainingError = FALSE){
  methodstitle = sapply(methods, function(cl) cl$classname)
  p = length(methods)
  results = list()
  # results.selected10 = list() # 10% pvalue cutoff
  results.selected5 = list() # 5% pval cutoff
  results.selected1 = list() # 1% pval cutoff
  #selectRate10 = rep(0,B)
  selectRate5 = rep(0,B)
  selectRate1 = rep(0,B)
  for(i in 1:B){
   # if (i%%20 == 0) print(i) # show progress
    traindata = simulateData(trainSize, priorGen, sampleSizeGen, scale)
    if(trainingError) {
      testdata = traindata
    } else {
      testdata = simulateData(testSize, priorGen, sampleSizeGen, scale)
    }
    testdata = testdata %>% mutate(tstat = y_scaled*sqrt(sample_size),pval = 2*pnorm(abs(tstat),lower=FALSE))
    # selectRate10[i] = mean(testdata$pval<0.1)
    selectRate5[i] = mean(testdata$pval<0.05)
    selectRate1[i] = mean(testdata$pval<0.01)
    
    run = function(cl, traindata, testdata){
      rmse = rep(NA, 3)
      cicov = rep(NA, 3)
      vrrate = rep(NA, 3)
      mdl = cl$new(traindata)
      if(!is.null(mdl$includeVar)) {mdl$includeVar = TRUE}
      trainedModel = mdl$train()
      preds = trainedModel$predict(testdata)
      if(!is.atomic(preds)){
        fits = preds$fit
        if(!is.null(cl$asymCI) && cl$asymCI){
          ciupper = fits$ciupper
          cilower = fits$cilower
        }else{
          sdfits = sqrt(preds$varfit)
          # if(any(is.nan(sdfits))){
          #   print(preds$varfit)
          # }
          ciupper = fits+1.96*sdfits
          cilower = fits-1.96*sdfits
          varfit = preds$varfit
          #print(ciupper)
        }
      } else{
        fits = preds
        varfit = testdata$scale^2/testdata$sample_size
        ciupper = fits+1.96*sqrt(testdata$scale^2/testdata$sample_size)
        cilower = fits-1.96*sqrt(testdata$scale^2/testdata$sample_size)
      }
      rmse[1] = caret::RMSE(fits, testdata$mu)
      cicov[1] = mean(testdata$mu<= ciupper & testdata$mu >= cilower, na.rm=TRUE)
      vrrate[1] = mean(varfit/(testdata$scale^2/testdata$sample_size), na.rm=TRUE)
      
      selected = testdata$pval < 0.05
      rmse[2] = caret::RMSE(fits[selected], testdata$mu[selected])
      cicov[2] = mean((testdata$mu<= ciupper & testdata$mu >= cilower)[selected], na.rm=TRUE)
      vrrate[2] = mean(varfit[selected]/(testdata$scale^2/testdata$sample_size)[selected], na.rm=TRUE)
      
      selected = testdata$pval < 0.01
      rmse[3] = caret::RMSE(fits[selected], testdata$mu[selected])
      cicov[3] = mean((testdata$mu<= ciupper & testdata$mu >= cilower)[selected], na.rm=TRUE)
      vrrate[3] = mean(varfit[selected]/(testdata$scale^2/testdata$sample_size)[selected], na.rm=TRUE)
      
      #rmse2 = caret::RMSE(preds, testdata$mu) # another score if needed
      c(rev(rmse),rev(cicov),rev(vrrate))
    }
    results[[i]] = vapply(methods, purrr::partial(run, traindata = traindata, testdata = testdata),FUN.VALUE = rep(0,9))
  }
  report = apply(vapply(results, identity, FUN.VALUE = results[[1]]), c(1,2), mean, na.rm=TRUE)

  rownames(report) = c("p<0.01","p<0.05", "All", "p<0.01-CI","p<0.05-CI", "All-CI", "p<0.01-VR","p<0.05-VR", "All-VR")
  colnames(report) = methodstitle
  report = as_tibble(report,rownames = "Selection")

  report$selectionRate = rep(c(1, mean(selectRate5), mean(selectRate1)),3)

  report
}



