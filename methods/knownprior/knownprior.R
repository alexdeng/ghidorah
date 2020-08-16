library(R6)

# When EXACT prior is known. Implement normal, Huber, laplace and t prior when all parameters are KNOWN. 
# Here the R6 class generator is itself output of a function that takes parameters as arguments. 
KnownNormalPrior = function(sigma){
  R6Class("KnownNormalPrior", 
    public = list(
      includeVar = TRUE,
      initialize= function(trainingSet) {}, # no action needed
      train = function(){invisible(self)}, 
      predict = function(newdata){
        required = c("y_scaled","scale","sample_size")
        assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
        B = sigma^2/(sigma^2+1/newdata$sample_size)
        preds =  B * newdata$y_scaled*newdata$scale
        if(self$includeVar){
          varfit = B / newdata$sample_size * newdata$scale^2
          return(list(fit = preds, varfit = varfit))
        }
        return(preds)
      }
    )
    )
}




KnownLaplacePrior = function(nu){ #nu: standard deviation of laplace prior
  R6Class("KnownLaplacePrior", 
    public = list(
      includeVar = TRUE,
      initialize= function(trainingSet) {}, # no action needed
      train = function(){invisible(self)}, 
      predict = function(newdata){
        required = c("y_scaled","scale","sample_size")
        assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
        b = sqrt(2)/nu/newdata$sample_size
        cy = sqrt(2)/nu*newdata$y_scaled
        Fy = exp(cy)*pnorm(sqrt(newdata$sample_size)*(-newdata$y_scaled-b))
        Gy = exp(-cy)*pnorm(sqrt(newdata$sample_size)*(newdata$y_scaled-b))
        wy = Fy/(Fy+Gy)
        preds = 2*wy*b + newdata$y_scaled-b
        if(self$includeVar){
          fy = nu*sqrt(newdata$sample_size)/sqrt(2)*exp(cy)*dnorm(sqrt(newdata$sample_size)*(-newdata$y_scaled-b))
          laplaceVar = 1/newdata$sample_size - 4/newdata$sample_size^2/nu^2*((Fy+Gy)*fy-2*Fy*Gy)/(Fy+Gy)^2
          return(list(fit = preds*newdata$scale, varfit = laplaceVar * newdata$scale^2))
        }
        return(preds*newdata$scale)
      }
    )
  )
}


KnownHuberPrior = function(b, K) {
  R6Class("KnownHuberPriorEB",
          public = list(
            includeVar = TRUE,
            initialize= function(trainingSet){}, #no action needed   
            train = function(){invisible(self)},
            predict = function(newdata){
              #check whether input data has all needed columns
              required = c("y_scaled","scale","sample_size")
              assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
              
              #cimpute posterior mean
              testdata = newdata %>% (private$transformFeature)
              if(self$includeVar){
                return(list(fit = testdata$huberfit*testdata$scale, varfit = testdata$huberVar*testdata$scale^2))
              }
              return(testdata$huberfit*testdata$scale)
            }
          ),
          private = list(
           #function to do the actual heavylifting
            #replacing the "EBHuber" function in the original code
            transformFeature = function(df) {
              #transform data
              y = df$y_scaled
              sigmas = 1 / df$sample_size
              #get paramaters in place
              mu = 0
              
              #compute posterior mean and variance
              v = K/b
              mu.t = (y/sigmas + mu/v)/(1/sigmas + 1/v)
              mu.b = y - b*sigmas
              mu.bb = y + b*sigmas
              sigma.t = 1/(1/v + 1/sigmas)
              W = K*b/2*(1+sigmas/v)
              Phi1 =  pnorm((mu + K-mu.t)/sqrt(sigma.t))-pnorm((mu-K-mu.t)/sqrt(sigma.t))
              Phi2 = pnorm(-(mu + K-mu.b)/sqrt(sigmas))
              Phi3 = pnorm((mu - K - mu.bb)/sqrt(sigmas))
              phi1 = dnorm((mu + K - mu.t)/sqrt(sigma.t))
              phi2 =  dnorm((mu -K - mu.t)/sqrt(sigma.t))
              phi3 =  dnorm((mu +K - mu.b)/sqrt(sigmas))
              phi4 = dnorm((mu-K - mu.bb)/sqrt(sigmas))
              U = -(y - mu)^2/(2*(sigmas + v))    
              
              #calculate log(m(y))
              logI = 1/2*(log(sigma.t)-log(sigmas))  - (y-mu)^2/(2*(v+sigmas)) + log(Phi1)
              logII = b*(mu-y) + W  + log(Phi2)
              logIII = b*(-mu+y) + W + log(Phi3)
              logmax = pmax(logI, logII, logIII)
              logm = logmax + log(exp(logI-logmax) + exp(logII- logmax) + exp(logIII - logmax))
              
              #calculate m'/m
              I= sqrt(sigma.t/sigmas)*exp(U)*(sqrt(sigma.t)/sigmas*(phi2-phi1) - (y-mu)/(sigmas + v)*Phi1)
              m1.p = which(phi3/sqrt(sigmas) > b*Phi2)
              m1.n = which(phi3/sqrt(sigmas) <= b*Phi2)
              m2.p = which(phi4/sqrt(sigmas) < b*Phi3)
              m2.n = which(phi4/sqrt(sigmas) >= b*Phi3)
              
              II = c()
              II[m1.p] = W[m1.p] + b*(mu-y[m1.p]) + log(phi3[m1.p]/sqrt(sigmas[m1.p]) - b*Phi2[m1.p])
              II[m1.n] = W[m1.n] + b*(mu-y[m1.n]) + log(-phi3[m1.n]/sqrt(sigmas[m1.n]) + b*Phi2[m1.n])
              
              III = c()
              III[m2.p] = W[m2.p] + b*(y[m2.p] - mu) + log(b*Phi3[m2.p] - phi4[m2.p]/sqrt(sigmas[m2.p]))
              III[m2.n] = W[m2.n] + b*(y[m2.n] - mu) + log(phi4[m2.n]/sqrt(sigmas[m2.n]) - b*Phi3[m2.n])
              
              P1 = I/exp(logm) + ((phi3/sqrt(sigmas) > b*Phi2)*exp(II - logm) - (phi3/sqrt(sigmas) <= b*Phi2)*(exp(II - logm))) + ((phi4/sqrt(sigmas) <b*Phi3)*exp(III - logm) - (phi4/sqrt(sigmas) >= b*Phi3)*(exp(III - logm)))
              
              #calculate m''/m
              C1 = ((y-mu)^2 - (sigmas + v))/(sigmas + v)^2
              I = sqrt(sigma.t/sigmas)*exp(U)*( C1*Phi1 + sqrt(sigma.t)/sigmas*(2*(y-mu)/(sigmas + v)-(mu+K-mu.t)/sigmas)*phi1 +
                                                  sqrt(sigma.t)/sigmas*(-2*(y-mu)/(sigmas + v)+(mu-K-mu.t)/sigmas )*phi2)
              C2 = ((mu+K-mu.b)/sigmas - 2*b)/sqrt(sigmas)
              C3 = (-(mu-K-mu.bb)/sigmas - 2*b)/sqrt(sigmas)
              m1.p = which(C2*phi3 + b^2*Phi2>=0)
              m1.n = which(C2*phi3+ b^2*Phi2< 0)
              m2.p = which(C3*phi4+ b^2*Phi3>=0)
              m2.n = which(C3*phi4+ b^2*Phi3<0)
              
              II = c()
              II[m1.p] = W[m1.p] + b*(mu-y[m1.p]) + log(C2[m1.p]*phi3[m1.p] +b^2*Phi2[m1.p])
              II[m1.n] = W[m1.n] + b*(mu-y[m1.n]) + log(-C2[m1.n]*phi3[m1.n] -b^2*Phi2[m1.n])
              
              III = c()
              III[m2.p] = W[m2.p] + b*(y[m2.p] - mu) + log(C3[m2.p]*phi4[m2.p] + b^2*Phi3[m2.p])
              III[m2.n] = W[m2.n] + b*(y[m2.n] - mu) + log(-C3[m2.n]*phi4[m2.n] - b^2*Phi3[m2.n])
              
              P2 =I/exp(logm) + ((C2*phi3+ b^2*Phi2>=0)*exp(II - logm) - (C2*phi3+ b^2*Phi2<0)*(exp(II - logm))) + ((C3*phi4+ b^2*Phi3>=0)*exp(III - logm) - (C3*phi4+ b^2*Phi3<0)*(exp(III - logm)))
              
              #final output
              df %>% mutate(#posterior mean
                huberfit = y_scaled + sigmas * P1,
                huberfit = huberfit*(sign(y_scaled)==sign(huberfit)), 
                #posterior variance
                huberVar = sigmas + sigmas^2*(P2 - P1^2),
                huberVar = if_else(is.nan(huberVar), 1/sample_size, huberVar))
            }
          )
  )
}


#' p is the proportion of zeros
KnownZeroInflatedNormalPrior = function(sigma, p = 0){
  assertthat::assert_that(p<1, msg = "proportion of 0 cannot be 1!")
  R6Class("KnownZeroInflatedNormalPrior", 
    public = list(
      includeVar = TRUE,
      initialize= function(trainingSet) {}, # no action needed
      train = function(){invisible(self)}, 
      predict = function(newdata){
        required = c("y_scaled","scale","sample_size")
        assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
        preds = sigma^2/(sigma^2+1/newdata$sample_size) * newdata$y_scaled
        postvar = sigma^2/(sigma^2+1/newdata$sample_size) / newdata$sample_size * newdata$scale^2
        lr = dnorm(newdata$y_scaled, 0, sqrt(sigma^2+1/newdata$sample_size))/dnorm(newdata$y_scaled, 0, sqrt(1/newdata$sample_size))
        postodds = (1-p)/p*lr
        posth1 = postodds/(1+postodds)
        fit = preds*newdata$scale*posth1
        varfit = posth1 * postvar + posth1*(1-posth1)*fit^2
        if(self$includeVar){
          return(list(fit = fit, varfit = varfit))
        }
        return(fit)
      }
    )
  )
}

#' p is the proportion of zeros
KnownZeroInflatedTPrior = function(alpha, tau, p = 0){
  assertthat::assert_that(p<1, msg = "proportion of 0 cannot be 1!")
  R6Class("KnownZeroInflatedTPrior", 
          public = list(
            includeVar = TRUE,
            initialize= function(trainingSet) {}, # no action needed
            train = function(){invisible(self)}, 
            predict = function(newdata){
              required = c("y_scaled","scale","sample_size")
              assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
              adjs = 1/newdata$sample_size*(alpha+1)*newdata$y_scaled / (alpha*tau^2+newdata$y_scaled^2)
              signs = sign(newdata$y_scaled)
              preds = newdata$y_scaled - adjs
              overadjed = (preds*newdata$y_scaled<0)
              preds[overadjed] = 0
              lr = private$my(newdata$y_scaled,newdata$sample_size)/dnorm(newdata$y_scaled, 0, sqrt(1/newdata$sample_size))
              lr[lr<0] = 0 #numeric issue in approximation
              lr[lr>10000] = 10000 #bound by 1000
              postodds = (1-p)/p*lr
              posth1 = postodds/(1+postodds)
              fit = posth1*preds*newdata$scale
              if(self$includeVar){
                #varadj = private$varadj(newdata$y_scaled, newdata$sample_size)
                #var_scaled = 1/newdata$sample_size - varadj
                #var_scaled[var_scaled<0] = 1/newdata$sample_size[var_scaled<0]
                #varfit = posth1*newdata$scale^2*var_scaled + posth1*(1-posth1)*fit^2  #t posterior variance is set unshrinked
                # if(any(varfit<0)){
                #   ind = varfit<0
                #   print(cbind(var_scaled[ind],posth1[ind],fit[ind],1/newdata$sample_size[ind]))
                # }
                return(list(fit = fit,varfit = 1/newdata$sample_size*newdata$scale^2 )) # varfit todo
              }
              return(fit)
            }
            # varadj2 =  function(y,n){
            #   ty0 = private$ty(y,0)
            #   ty1 = private$ty(y,1)
            #   ty2 = private$ty(y,2)
            #   (y^2*ty1^2/tau^4/ty0^2-(y^2*(alpha+2)*ty2/tau^4/alpha-ty1/tau^2)/ty0)/n^2
            # }
          ),
          private = list(
            ty = function(y,m){
              dt(y/sqrt(tau^2*alpha/(alpha+2*m)),df=alpha+2*m,ncp=0)
            },
            varadj = function(y,n){
              ty0 = private$ty(y,0)
              ty1 = private$ty(y,1)
              ty2 = private$ty(y,2)
              (y^2*ty1^2/tau^4/ty0^2-(y^2*(alpha+2)*ty2/tau^4/alpha-ty1/tau^2)/ty0)/n^2
            },
            my = function(y,n){
              ty0 = private$ty(y,0)
              ty1 = private$ty(y,1)
              ty2 = private$ty(y,2)
              ty0 - 1/(2*n*tau^2)*ty1 + (y^2*(alpha+2))*ty2/(2*n*tau^4*alpha)
            }
          )
  )
}



KnownTPrior = function(alpha, tau){ #alpha: degree of freedom, tau: scale
  R6Class("KnownTPrior", 
          public = list(
            includeVar = TRUE,
            initialize= function(trainingSet) {}, # no action needed
            train = function(){invisible(self)} ,
            predict = function(newdata){
              required = c("y_scaled","scale","sample_size")
              assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
              adjs = 1/newdata$sample_size*(alpha+1)*newdata$y_scaled / (alpha*tau^2+newdata$y_scaled^2)
              signs = sign(newdata$y_scaled)
              preds = newdata$y_scaled - adjs
              overadjed = (preds*newdata$y_scaled<0)
              preds[overadjed] = 0
              if(self$includeVar){
                varadj = private$varadj(newdata$y_scaled, newdata$sample_size)
                var_scaled = 1/newdata$sample_size - varadj
                var_scaled[var_scaled<0] = 1/newdata$sample_size[var_scaled<0]*0.1 # bounded away from 0
                return(list(fit = preds*newdata$scale,varfit = var_scaled*newdata$scale^2 )) 
              }
              return(preds*newdata$scale)
             }
            # varadj2 =  function(y,n){
            #   ty0 = private$ty(y,0)
            #   ty1 = private$ty(y,1)
            #   ty2 = private$ty(y,2)
            #   (y^2*ty1^2/tau^4/ty0^2-(y^2*(alpha+2)*ty2/tau^4/alpha-ty1/tau^2)/ty0)/n^2
            # },
            # meanadj = function(y,n){
            #   1/n*(alpha+1)*y / (alpha*tau^2+y^2)
            # }
          ),
          private = list(
            ty = function(y,m){
              dt(y/sqrt(tau^2*alpha/(alpha+2*m)),df=alpha+2*m,ncp=0)
            },
            varadj = function(y,n){
              ty0 = private$ty(y,0)
              ty1 = private$ty(y,1)
              ty2 = private$ty(y,2)
              (y^2*ty1^2/tau^4/ty0^2-(y^2*(alpha+2)*ty2/tau^4/alpha-ty1/tau^2)/ty0)/n^2
            }
          )
  )
}
