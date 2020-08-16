library(R6)
source("methods/MixturePrior/mixture.R")

MACRO_USEFQVAR = TRUE

# this is linear model using predictors inspired by the true posterior form for t-prior (df 3,5,10,20), laplace prior and normal prior(J-S). 
# hope is this will cover a broad family of priors
TARwESSplitting = R6Class("TARwES",
  public = list(
    nu = NULL, # laplace sd
    includeVar = FALSE,
    USEFQVAR = MACRO_USEFQVAR,
    SUREnu = TRUE, # Use SURE to fit standard deviation of prior -- nu
    initialize= function(trainingSet){
      required = c("y_scaled", "y_scaled_A","y_scaled_B", "sample_size_A", "sample_size")
      assertthat::assert_that(length(intersect(required,colnames(trainingSet))) == length(required), msg = "Input does not have all the required columns: y_scaled, y_scaled_A, y_scaled_B, sample_size_A, sample_size!")
      # the estimated scale parameter nu for laplace prior
      total.var = var(trainingSet$y_scaled)  
      noise.var = mean(1/trainingSet$sample_size)
      prior.var = max(total.var-noise.var, noise.var * 0.1) # if numerical value is negative, bound by 10% of noise's variance
      if(self$SUREnu){
        surerisk = function(nu){
          B = 1/trainingSet$sample_size/(1/trainingSet$sample_size + nu^2)
          mean(B^2*trainingSet$y_scaled^2+2*nu^2*B)
        }
        self$nu = max(optim(sqrt(prior.var),surerisk, method="L-BFGS-B", lower=0, upper=sqrt(total.var))$par,sqrt(0.1*noise.var)) #bounded by 10% of noise var
      }else{
        self$nu = sqrt(prior.var)
      }
      private$training = trainingSet %>% select(required) %>% select(-y_scaled)
    },
    train = function(){
      enriched = (private$training) %>% (private$enrichTrainingData) 
      modx = enriched %>% select(js, lapfit) %>% as.matrix
      W = diag(sqrt(enriched$sample_size/1000))
      private$model = nnls::nnls(W%*%modx, W%*%enriched$y_scaled_O) # weighted nonnegtive regression for regularization (L2 regularization gave similar performance)
      if(self$includeVar){
        resids = as.vector(residuals(private$model))
        varfeatures = enriched %>% (private$varFeature) %>% mutate(residsqAdj = resids^2- 1/sample_size)
        private$varmodel = lm(residsqAdj ~ -1 + fqvar  + jsVar + laplaceVar, data = varfeatures, weights = enriched$sample_size/1000) 
      }
      invisible(self)
    },
    predict = function(newdata){
      required = c("y_scaled","scale","sample_size")
      assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
      testdata = newdata %>% (private$enrichTestingData)
      newx = testdata %>% select(js, lapfit) %>% as.matrix
      scaledPred = as.vector(newx %*% (coef(private$model)))
      pred = scaledPred * newdata$scale * (sign(scaledPred)==sign(testdata$y_scaled))
      if(self$includeVar){
        vardata = testdata %>% (private$varFeature)
        if(self$USEFQVAR){
          varfit = vardata$fqvar
        }else{
          scaledvarPred = predict(private$varmodel, newdata = vardata) 
          regfit = scaledvarPred*newdata$scale^2
          js = newdata$scale^2*vardata$jsVar
          #laplace = newdata$scale^2*vardata$laplaceVar
          varfit = ifelse(regfit/js < 0.1, js, regfit)
        }
          return(list(fit = pred, varfit = varfit, jsfit = testdata$js*newdata$scale, lapfit = testdata$lapfit*newdata$scale))
        }
      return(pred)
    },
    varEst = function(newdata){
      required = c("y_scaled","scale","sample_size")
      assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
      testdata = newdata %>% (private$enrichTestingData)
      scaledPred = predict(private$model, newdata = testdata)
      pred = scaledPred * newdata$scale * (sign(scaledPred)==sign(testdata$y_scaled))
      vardata = testdata %>% (private$varFeature)
      regpred = newdata$scale^2*(predict(private$varmodel, newdata = vardata))
      js = newdata$scale^2*vardata$jsVar
      laplace = newdata$scale^2*vardata$laplaceVar
      return(tibble(reg = regpred, js = js, laplace = laplace, fq = 1/newdata$sample_size))
    },
    printModel = function() print(private$model),
    printVarModel = function() print(private$varmodel)
  ),                        
  private = list(
    training = NULL,
    model = NULL,
    varmodel = NULL,
    enrichTrainingData = function(df){
      df %>% 
        rename(y_scaled=y_scaled_A, y_scaled_O=y_scaled_B) %>% mutate(sample_size = sample_size_A) %>% (private$symmetrify) %>%  # sample size for a split is a half
        (private$transformFeature)
    },
    symmetrify = function(df){
      bind_rows(df, df %>% mutate(y_scaled_O = - y_scaled_O, y_scaled = - y_scaled))
    },
    enrichTestingData = function(df){
      df %>% (private$transformFeature) 
    },
    transformFeature = function(df) {
      nu = self$nu
      assertthat::assert_that(!is.null(nu), msg="Parameter nu has not been initiated from training data!")
      # sigma is 1 after scaled
      tmp = df %>% mutate(b = 1/sample_size/nu*sqrt(2), cy = sqrt(2)/nu*(y_scaled), Fy = exp(cy)*pnorm(sqrt(sample_size)*(-y_scaled-b)), Gy = exp(-cy)*pnorm(sqrt(sample_size)*(y_scaled-b)), wy = Fy/(Fy+Gy), wy = if_else(is.nan(wy),(1-sign(y_scaled))/2,wy), wyb = wy*b, js = nu^2/(1/sample_size+nu^2)*y_scaled, lapfit = y_scaled + 2*wyb - b, lapfit = lapfit*(sign(y_scaled)==sign(lapfit)))
      #tmp %>% mutate(tf3 = private$tfeature(y_scaled, 3, nu)/sample_size, tf5 = private$tfeature(y_scaled, 5, nu)/sample_size, tf10 = private$tfeature(y_scaled,10, nu)/sample_size, tf20 = private$tfeature(y_scaled, 20, nu)/sample_size)
    },
    tfeature = function(y,degfree,nu) {
      tausq = nu^2 *(degfree-1)/degfree
      (degfree+1)*y/(degfree*tausq+y^2)
    }, 
    varFeature = function(df){
      nu = self$nu
      df %>% mutate(fqvar = 1/sample_size,jsVar = nu^2/(1/sample_size+nu^2)/sample_size,fy = nu*sqrt(sample_size)/sqrt(2)*exp(cy)*dnorm(sqrt(sample_size)*(-y_scaled-b)), laplaceVar = 1/sample_size - 4/sample_size^2/nu^2*((Fy+Gy)*fy-2*Fy*Gy)/(Fy+Gy)^2, laplaceVar = if_else(is.nan(laplaceVar),1/sample_size, laplaceVar))
    }
  ),
  active = list(
    trained = function(){
      !is.null(private$model)
    },
    postMeanModel = function(){private$model},
    postVarModel = function(){private$varmodel},
    trainingData = function(){private$training %>% (private$enrichTrainingData)  }
  )
)


# this is linear model using predictors inspired by the true posterior form for t-prior (df 3,5,10,20), laplace prior and normal prior(J-S). 
# hope is this will cover a broad family of priors
GhidorahSplitting = R6Class("TARwES+",
  public = list(
    nu = NULL, # laplace sd
    includeVar = FALSE,
    USEFQVAR = MACRO_USEFQVAR,
    Ghidorah = NULL,
    SUREnu = TRUE, # Use SURE to fit standard deviation of prior -- nu
    initialize= function(trainingSet){
      required = c("y_scaled", "y_scaled_A","y_scaled_B", "sample_size_A", "sample_size","scale")
      assertthat::assert_that(length(intersect(required,colnames(trainingSet))) == length(required), msg = "Input does not have all the required columns: y_scaled, y_scaled_A, y_scaled_B, sample_size_A, sample_size!")
      # initialize the inner Ghidorah model
      self$Ghidorah = ZeroNormalLaplaceMixEB$new(trainingSet)
      self$Ghidorah$includeVar = FALSE
      self$Ghidorah$train()
      # the estimated scale parameter nu for laplace prior
      total.var = var(trainingSet$y_scaled)  
      noise.var = mean(1/trainingSet$sample_size)
      prior.var = max(total.var-noise.var, noise.var * 0.1) # if numerical value is negative, bound by 10% of noise's variance
      if(self$SUREnu){
        surerisk = function(nu){
          B = 1/trainingSet$sample_size/(1/trainingSet$sample_size + nu^2)
          mean(B^2*trainingSet$y_scaled^2+2*nu^2*B)
        }
        self$nu = max(optim(sqrt(prior.var),surerisk, method="L-BFGS-B", lower=0, upper=sqrt(total.var))$par,sqrt(0.1*noise.var)) #bounded by 10% of noise var
      }else{
        self$nu = sqrt(prior.var)
      }
      private$training = trainingSet %>% select(required) %>% select(-y_scaled)
    },
    train = function(){
      enriched = (private$training) %>% (private$enrichTrainingData) 
      gfit = self$Ghidorah$predict(enriched)/enriched$scale
      modx = enriched %>%  mutate(gfit = gfit) %>% select(y_scaled,js, lapfit, gfit) %>% as.matrix
      W = diag(sqrt(enriched$sample_size/1000))
      private$model = nnls::nnls(W%*%modx, W%*%enriched$y_scaled_O) # weighted nonnegtive regression. L2 regularization gave similar performance
      if(self$includeVar){
        resids = as.vector(residuals(private$model))
        varfeatures = enriched %>% (private$varFeature) %>% mutate(residsqAdj = resids^2- 1/sample_size)
        private$varmodel = lm(residsqAdj ~ -1 + fqvar  + jsVar + laplaceVar, data = varfeatures, weights = enriched$sample_size/1000) 
      }
      invisible(self)
    },
    predict = function(newdata){
      required = c("y_scaled","scale","sample_size")
      assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
      testdata = newdata %>% (private$enrichTestingData)
      gfit = self$Ghidorah$predict(newdata)/newdata$scale
      newx = testdata %>% mutate(gfit=gfit) %>% select(y_scaled,js, lapfit, gfit) %>% as.matrix
      scaledPred = as.vector(newx %*% (coef(private$model)))
      pred = scaledPred * newdata$scale * (sign(scaledPred)==sign(testdata$y_scaled))
      if(self$includeVar){
        vardata = testdata %>% (private$varFeature)
        if(self$USEFQVAR){
          varfit = vardata$fqvar
        }else{
          scaledvarPred = predict(private$varmodel, newdata = vardata) 
          regfit = scaledvarPred*newdata$scale^2
          js = newdata$scale^2*vardata$jsVar
          #laplace = newdata$scale^2*vardata$laplaceVar
          varfit = ifelse(regfit/js < 0.1, js, regfit)
        }
        return(list(fit = pred, varfit = varfit, jsfit = testdata$js*newdata$scale, lapfit = testdata$lapfit*newdata$scale))
      }
      return(pred)
    },
    varEst = function(newdata){
      required = c("y_scaled","scale","sample_size")
      assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
      testdata = newdata %>% (private$enrichTestingData)
      scaledPred = predict(private$model, newdata = testdata)
      pred = scaledPred * newdata$scale * (sign(scaledPred)==sign(testdata$y_scaled))
      vardata = testdata %>% (private$varFeature)
      regpred = newdata$scale^2*(predict(private$varmodel, newdata = vardata))
      js = newdata$scale^2*vardata$jsVar
      laplace = newdata$scale^2*vardata$laplaceVar
      return(tibble(reg = regpred, js = js, laplace = laplace, fq = 1/newdata$sample_size))
    },
    printModel = function() print(private$model),
    printVarModel = function() print(private$varmodel)
  ),                        
  private = list(
    training = NULL,
    model = NULL,
    varmodel = NULL,
    enrichTrainingData = function(df){
      df %>% 
        rename(y_scaled=y_scaled_A, y_scaled_O=y_scaled_B) %>% mutate(sample_size = sample_size_A) %>% (private$symmetrify) %>%  # sample size for a split is a half
        (private$transformFeature)
    },
    symmetrify = function(df){
      bind_rows(df, df %>% mutate(y_scaled_O = - y_scaled_O, y_scaled = - y_scaled))
    },
    enrichTestingData = function(df){
      df %>% (private$transformFeature) 
    },
    transformFeature = function(df) {
      nu = self$nu
      assertthat::assert_that(!is.null(nu), msg="Parameter nu has not been initiated from training data!")
      # sigma is 1 after scaled
      tmp = df %>% mutate(b = 1/sample_size/nu*sqrt(2), cy = sqrt(2)/nu*(y_scaled), Fy = exp(cy)*pnorm(sqrt(sample_size)*(-y_scaled-b)), Gy = exp(-cy)*pnorm(sqrt(sample_size)*(y_scaled-b)), wy = Fy/(Fy+Gy), wy = if_else(is.nan(wy),(1-sign(y_scaled))/2,wy), wyb = wy*b, js = nu^2/(1/sample_size+nu^2)*y_scaled, lapfit = y_scaled + 2*wyb - b, lapfit = lapfit*(sign(y_scaled)==sign(lapfit)))
      #tmp %>% mutate(tf3 = private$tfeature(y_scaled, 3, nu)/sample_size, tf5 = private$tfeature(y_scaled, 5, nu)/sample_size, tf10 = private$tfeature(y_scaled,10, nu)/sample_size, tf20 = private$tfeature(y_scaled, 20, nu)/sample_size)
    },
    tfeature = function(y,degfree,nu) {
      tausq = nu^2 *(degfree-1)/degfree
      (degfree+1)*y/(degfree*tausq+y^2)
    }, 
    varFeature = function(df){
      nu = self$nu
      df %>% mutate(fqvar = 1/sample_size,jsVar = nu^2/(1/sample_size+nu^2)/sample_size,fy = nu*sqrt(sample_size)/sqrt(2)*exp(cy)*dnorm(sqrt(sample_size)*(-y_scaled-b)), laplaceVar = 1/sample_size - 4/sample_size^2/nu^2*((Fy+Gy)*fy-2*Fy*Gy)/(Fy+Gy)^2, laplaceVar = if_else(is.nan(laplaceVar),1/sample_size, laplaceVar))
    }
  ),
  active = list(
    trained = function(){
      !is.null(private$model)
    },
    postMeanModel = function(){private$model},
    postVarModel = function(){private$varmodel},
    trainingData = function(){private$training %>% (private$enrichTrainingData)  }
  )
)




