HuberPriorEB = 
  R6Class("HuberPriorEB",
          public = list(
            #Huber distribution parameters
            mu = 0, #always symmetrize the data so mu = 0 automatically
            b = NULL,
            K = NULL,
            includeVar = FALSE,
            #use SURE to fit parameters of Huber prior
            SURE = FALSE,
            initialize= function(trainingSet){
              #check whether input data has all needed columns
              required = c("y_scaled","scale", "sample_size")
              assertthat::assert_that(length(intersect(required,colnames(trainingSet))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
              #using SURE method
              if(self$SURE) {
                para.huber.opt = optim(par = c(var(trainingSet$y_scaled), 1), lower = c(1e-4, 1e-4), upper = c(1e6, 1e6), method="L-BFGS-B", gr = NULL, private$Huber.SURE, df = trainingSet)$par
                self$b = para.huber.opt[1]
                self$K = para.huber.opt[2]
              }
              #using marginal MLE method
              else {
                para.huber.opt = optim(par = c(var(trainingSet$y_scaled), 1), lower = c(1e-4, 1e-4), upper = c(1e6, 1e6), method="L-BFGS-B", gr = NULL, private$Huber.negloglik, df = trainingSet)$par
                self$b = para.huber.opt[1]
                self$K = para.huber.opt[2]
              }
            },   
            train = function(){invisible(self)},
            predict = function(newdata){
              #check whether input data has all needed columns
              required = c("y_scaled","scale","sample_size")
              assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
              
              #check whether parameters of prior have been learned
              assertthat::assert_that(!is.null(self$mu), msg="Parameter mu has not been initiated from training data!")
              assertthat::assert_that(!is.null(self$b), msg="Parameter b has not been initiated from training data!")
              assertthat::assert_that(!is.null(self$K), msg="Parameter K has not been initiated from training data!")
              
              #cimpute posterior mean and variance (optional)
              testdata = newdata %>% (private$transformFeature)
              if(self$includeVar){
                return(list(fit = testdata$huberfit*testdata$scale, varfit = testdata$huberVar*testdata$scale^2))
              }
              return(testdata$huberfit*testdata$scale)
            },
            debugdata = function(newdata){
              newdata %>% (private$transformFeature)
            }
          ),
          private = list(
            #function to compute SURE
            Huber.SURE = function(para, df){
              #transform data
              y = df$y_scaled
              sigmas = 1 / df$sample_size
              
              #get parameters in place
              mu = 0 #always symmetrize the data so mu = 0 automatically
              b = para[1]
              K = para[2]
              
              v = K/b
              mu.t = (y/sigmas + mu/v)/(1/sigmas + 1/v)
              mu.b = y - b*sigmas
              mu.bb = y + b*sigmas
              sigma.t = 1/(1/v + 1/sigmas)
              W = K*b/2*(1+sigmas/v)
              Phi1 =  pnorm((mu+K-mu.t)/sqrt(sigma.t))-pnorm((mu-K-mu.t)/sqrt(sigma.t))
              Phi2 = pnorm(-(mu+K-mu.b)/sqrt(sigmas))
              Phi3 = pnorm((mu-K-mu.bb)/sqrt(sigmas))
              phi1 = dnorm((mu +K - mu.t)/sqrt(sigma.t))
              phi2 =  dnorm((mu -K - mu.t)/sqrt(sigma.t))
              phi3 =  dnorm((mu +K - mu.b)/sqrt(sigmas))
              phi4 = dnorm((mu-K - mu.bb)/sqrt(sigmas))
              U = -(y-mu)^2/(2*(sigmas+ v))    
              
              
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
              
              
              #put together
              sure = sum(sigmas^2*(2*P2 - P1^2))
              return(sure)
            },
            
            #function to compute negative log joint marginal likelihood under Huber prior
            Huber.negloglik = function(paras, df){
              #transform data
              y = df$y_scaled
              sigmas = 1 / df$sample_size
              
              #get parameters in place
              mu = 0 #always symmetrize the data so mu = 0 automatically
              b = paras[1]
              K = paras[2]
              
              mu.t = (y/sigmas + mu*b/K)/(1/sigmas + b/K)
              mu.b = y - b*sigmas
              mu.bb = y + b*sigmas
              sigma.t = 1/(b/K + 1/sigmas)
              W = b*K*(1+b*sigmas/K)/2
              phi1 =  pnorm((mu+K-mu.t)/sqrt(sigma.t)) - pnorm((mu-K-mu.t)/sqrt(sigma.t))
              logphi2 = pnorm(-(mu+K-mu.b)/sqrt(sigmas), log.p = T)
              logphi3 = pnorm((mu-K-mu.bb)/sqrt(sigmas), log.p = T)
              C = sqrt(2*pi*K/b)*(pnorm(K/sqrt(K/b)) - pnorm(-K/sqrt(K/b))) + 2*exp(-K*b/2)/b
              logI = 1/2*(log(sigma.t)-log(sigmas))  -(y-mu)^2/(2*(K/b+sigmas)) + log(phi1)
              logII = b*(mu-y) + W + logphi2
              logIII = b*(-mu+y) +W + logphi3
              logmax = pmax(logI, logII, logIII)
              neglog = -sum(logmax+log(exp(logI - logmax) + exp(logII - logmax) + exp(logIII - logmax))) + length(y)*log(C)
              return(neglog)
            },
            
            #function to do the actual heavylifting
            #replacing the "EBHuber" function in the original code
            transformFeature = function(df) {
              #transform data
              y = df$y_scaled
              sigmas = 1 / df$sample_size
              #get paramaters in place
              mu = self$mu
              b = self$b
              K = self$K
              #check whether parameters of prior have been learned
              assertthat::assert_that(!is.null(self$mu), msg="Parameter mu has not been initiated from training data!")
              assertthat::assert_that(!is.null(self$b), msg="Parameter b has not been initiated from training data!")
              assertthat::assert_that(!is.null(self$K), msg="Parameter K has not been initiated from training data!")
              
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
