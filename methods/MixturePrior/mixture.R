library(R6)
# When this is set to be TRUE, Ghidorah will cap variance by MLE's variance sigma^2/n, i.e. variance will only reduce. 
MACRO_CAPVAR = TRUE


# three group, zero, normal and laplace
ZeroNormalLaplaceMixEB = 
  R6Class("Ghidorah",
          public = list(
            nu = NULL, # a vector of sd for normal and laplace prior
            p = NULL, # a vector of prior probability for zero, normal and laplace
            SUREnu = TRUE, # Use SURE to fit standard deviation of prior -- nu
            includeVar = TRUE,
            capVar = MACRO_CAPVAR,
            maxpZero = 0.99, # cap prior of zero to be 95% to avoid too overfit to zero
            maxIter = 5000, # max iteration for EM algorithm
            initialize= function(trainingSet){
              required = c("y_scaled","scale", "sample_size")
              assertthat::assert_that(length(intersect(required,colnames(trainingSet))) == length(required), msg = "Input does not have all the required columns: y_scaled_A, y_scaled_B, scale, sample_size!")
              # the estimated scale parameter nu for laplace prior
              total.var = mean(trainingSet$y_scaled^2)  
              noise.var = mean(1/trainingSet$sample_size)
              prior.var = max(total.var-noise.var, noise.var * 0.1) # if numerical value is negative, bound by 1% of noise's variance
              if(self$SUREnu){
                surerisk = function(nu){
                  B = 1/trainingSet$sample_size/(1/trainingSet$sample_size + nu^2)
                  mean(B^2*trainingSet$y_scaled^2+2*nu^2*B)
                }
                nu = max(optim(sqrt(prior.var),surerisk, method="L-BFGS-B", lower=0, upper=sqrt(total.var))$par,sqrt(0.1*noise.var)) #bounded by 10% of noise var
                private$nuCurr = c(nu, nu) #initial value
              }else{
                nu = sqrt(prior.var)
                private$nuCurr = c(nu, nu) #initial value
              }
              private$training = trainingSet %>% select(required) 
            },   
            # Use EM to fit the mixture prior model
            train = function(){
              enriched = (private$training) %>% (private$enrichTrainingData)
              for (i in 1:self$maxIter){
                private$EMIteration(enriched) # 10 iterations
                if(self$converged){
                  break
                }
              }
              invisible(self)
            },
            predict = function(newdata){
              required = c("y_scaled","scale","sample_size")
              assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
              assertthat::assert_that(!is.null(self$nu), msg="Parameter nu has not been trained from training data!")
              assertthat::assert_that(!is.null(self$p), msg="Parameter p has not been trained from training data!")
              testdata = newdata %>% (private$enrichTestingData)
              pred = (testdata$posteriorNormal*testdata$jsfit + testdata$posteriorLaplace*testdata$lapfit)*newdata$scale
              if(self$includeVar){
                vardata = testdata %>% (private$varFeature)
                postNorm = testdata$posteriorNormal
                postLap = testdata$posteriorLaplace
                fqvar = vardata$fqvar*newdata$scale^2
                varfit = (postNorm*vardata$jsVar + postLap*vardata$laplaceVar+postNorm*vardata$jsfit^2 + postLap*vardata$lapfit^2-(postNorm*testdata$jsfit + postLap*testdata$lapfit)^2)*newdata$scale^2
                if(self$capVar){
                  ind = varfit>fqvar #faster than ifelse
                  varfit[ind] = fqvar[ind] # cap at FQ variance. 
                }
                return(list(fit = pred, varfit = varfit, jsfit = vardata$jsfit*vardata$scale, lapfit = vardata$lapfit*vardata$scale, posteriorLaplace = postLap))
              }
              return(pred)
            },
            debugData = function(){
              (private$training) %>% (private$enrichTrainingData)
            },
            EMCurr = function(){
              list(p = private$pCurr, nu = private$nuCurr)
            }
          ),        
          private = list(
            pCurr = c(1/3,1/3,1/3), #current p vector (for iterative training)
            nuCurr = NULL, #current nu vector (for iterative training)
            training = NULL,
            enrichTrainingData = function(df){
              df %>% 
                (private$symmetrify)
            },
            #symmetrify the input data to force a symmetric prior
            symmetrify = function(df){
              bind_rows(df, df %>% mutate(y_scaled = - y_scaled))
            },
            enrichTestingData = function(df){
              df %>% (private$transformFeature) 
            },
            transformFeature = function(df) {
              assertthat::assert_that(!is.null(self$nu), msg="Parameter nu has not been initiated from training data!")
              assertthat::assert_that(!is.null(self$p), msg="Parameter p has not been initiated from training data!")
              normnu = self$nu[1]
              lapnu = self$nu[2]
              p = self$p
              # sigma is 1 after scaled
              df %>% mutate(a=exp(1/lapnu^2/sample_size) ,b = 1/sample_size/lapnu*sqrt(2), cy = sqrt(2)/lapnu*(y_scaled), Fy = exp(cy)*pnorm(sqrt(sample_size)*(-y_scaled-b)), Gy = exp(-cy)*pnorm(sqrt(sample_size)*(y_scaled-b)), wy = Fy/(Fy+Gy), wy = if_else(is.nan(wy),(1-sign(y_scaled))/2,wy), wyb = wy*b, jsfit_prime = normnu^2/(1/sample_size+normnu^2), jsfit = jsfit_prime*y_scaled, lapfit = y_scaled + 2*wyb - b, lapfit = lapfit*(sign(y_scaled)==sign(lapfit)), fy = lapnu*sqrt(sample_size)/sqrt(2)*exp(cy)*dnorm(sqrt(sample_size)*(-y_scaled-b)), w_prime = -((Fy+Gy)*fy-2*Fy*Gy)/(Fy+Gy)^2*sqrt(2)/lapnu, lapMarginal = a/lapnu/sqrt(2)*(Fy+Gy)*p[3], normMarginal = dnorm(y_scaled,0,sqrt(1/sample_size+normnu^2))*p[2], zeroMarginal = dnorm(y_scaled,0,sqrt(1/sample_size))*p[1], sumMarginal = zeroMarginal+normMarginal+lapMarginal, posteriorZero = zeroMarginal/sumMarginal, posteriorNormal = normMarginal/sumMarginal, posteriorLaplace = ifelse(is.nan(lapMarginal),1,lapMarginal/sumMarginal))
            },
            varFeature = function(df){
              lapnu = self$nu[2]
              df %>% mutate(fqvar = 1/sample_size,jsVar = jsfit_prime/sample_size,laplaceVar = 1/sample_size + 2*sqrt(2)/sample_size^2/lapnu*w_prime, laplaceVar = if_else(is.nan(laplaceVar),1/sample_size, laplaceVar))
            },
            MMVarEstimate = function(w,df){
              total.var = private$safeWeightedMean(df$y_scaled^2,w)  
              noise.var = private$safeWeightedMean(1/df$sample_size,w)
              prior.var = max(total.var-noise.var, noise.var * 0.2) # if numerical value is negative, bound by 20% of noise's variance
              sqrt(prior.var)
            },
            EMIteration = function(df){
              self$nu = private$nuCurr
              self$p = private$pCurr
              #update 
              logzeroMarginalLikelihood = dnorm(df$y_scaled,0,sqrt(1/df$sample_size),log = TRUE)+log(self$p[1])
              lognormalMarginalLikelihood = dnorm(df$y_scaled,0,sqrt(1/df$sample_size + self$nu[1]^2),log = TRUE)+log(self$p[2])
              zeroNormalRatio = exp(logzeroMarginalLikelihood - lognormalMarginalLikelihood)
              
              #lapace marginal likelihood is a/nu/sqrt(2)*(Fy+Gy) See Pericchi and Smith
              lapnu = self$nu[2]
              b = 1/df$sample_size/lapnu*sqrt(2)
              loga = 1/lapnu^2/df$sample_size
              cy = sqrt(2)/lapnu*(df$y_scaled)
              logFy = cy+pnorm(sqrt(df$sample_size)*(-df$y_scaled-b), log.p = TRUE)
              logGy = -cy+pnorm(sqrt(df$sample_size)*(df$y_scaled-b), log.p=TRUE)
              laplaceNormalRatio = 1/lapnu*(exp(logFy+loga - lognormalMarginalLikelihood)+exp(logGy+loga - lognormalMarginalLikelihood))/sqrt(2)*self$p[3]
              #print(list(p = self$p, nu = self$nu,lvn = laplaceNormalRatio[1:20], zvn = zeroNormalRatio[1:20]))

              #readline("press")
              laplaceNormalRatio[!is.finite(laplaceNormalRatio)] = 1e6
              zeroNormalRatio[!is.finite(zeroNormalRatio)] = 1e6
              tmp = 1 + zeroNormalRatio + laplaceNormalRatio
              posteriorZero = zeroNormalRatio/tmp
              posteriorNormal = 1/tmp
              posteriorLaplace = laplaceNormalRatio/tmp
              
              private$pCurr = c(mean(posteriorZero), mean(posteriorNormal), mean(posteriorLaplace))

              if(private$pCurr[1]>self$maxpZero) {
                private$pCurr = c(self$maxpZero, (1-self$maxpZero)/sum(private$pCurr[2:3])*private$pCurr[2:3])
              }
              private$nuCurr = c(private$MMVarEstimate(posteriorNormal, df), private$MMVarEstimate(posteriorLaplace,df))
            },
            safeWeightedMean = function(x,w){
              y = sum(w)
              if(y<.Machine$double.eps) return(mean(x))
              sum(x*w)/y
            }
          ),
          active = list(
            converged = function(){
              if(!is.null(self$p) & !is.null(self$nu)){
                dis1 = max(abs(self$p-private$pCurr))
                dis1<1e-3
              }else{
                FALSE
              }
            }
          )
  )

