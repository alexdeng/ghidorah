library(R6)

NormalPriorEB = 
  R6Class("NormalPriorEB(J-S)",
    public = list(
      nu = NULL, # laplace sd
      includeVar = FALSE,
      SUREnu = TRUE, # Use SURE to fit standard deviation of prior -- nu
      initialize= function(trainingSet){
        required = c("y_scaled","scale", "sample_size")
        assertthat::assert_that(length(intersect(required,colnames(trainingSet))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
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
      },   
      train = function(){invisible(self)},
      predict = function(newdata){
        required = c("y_scaled","scale","sample_size")
        assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
        assertthat::assert_that(!is.null(self$nu), msg="Parameter nu has not been initiated from training data!")
        shrinkage = self$nu^2/(self$nu^2+1/newdata$sample_size)
        if(self$includeVar){
          return(list(fit = shrinkage*newdata$y_scaled*newdata$scale, varfit = shrinkage*newdata$scale^2/newdata$sample_size))
        }
        return(shrinkage*newdata$y_scaled*newdata$scale)
      })
)


LaplacePriorEB = 
  R6Class("LaplacePriorEB",
    public = list(
      nu = NULL, # laplace sd
      includeVar = FALSE,
      SUREnu = TRUE, # Use SURE to fit standard deviation of prior -- nu
      initialize= function(trainingSet){
        required = c("y_scaled","scale", "sample_size")
        assertthat::assert_that(length(intersect(required,colnames(trainingSet))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
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
      },   
      train = function(){invisible(self)},
      predict = function(newdata){
        required = c("y_scaled","scale","sample_size")
        assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
        assertthat::assert_that(!is.null(self$nu), msg="Parameter nu has not been initiated from training data!")
        testdata = newdata %>% (private$transformFeature)
        if(self$includeVar){
          return(list(fit = testdata$lapfit*testdata$scale, varfit = testdata$laplaceVar*testdata$scale^2))
        }
        return(testdata$lapfit*testdata$scale)
      },
      debugdata = function(newdata){
        newdata %>% (private$transformFeature)
      }
      ),
    private = list(
      transformFeature = function(df) {
        nu = self$nu
        assertthat::assert_that(!is.null(nu), msg="Parameter nu has not been initiated from training data!")
        # sigma is 1 after scaled
        df %>% mutate(b = 1/sample_size/nu*sqrt(2), cy = sqrt(2)/nu*(y_scaled), Fy = exp(cy)*pnorm(sqrt(sample_size)*(-y_scaled-b)), Gy = exp(-cy)*pnorm(sqrt(sample_size)*(y_scaled-b)), wy = Fy/(Fy+Gy), wy = if_else(is.nan(wy),(1-sign(y_scaled))/2,wy), wyb = wy*b, lapfit = y_scaled + 2*wyb - b, lapfit = lapfit*(sign(y_scaled)==sign(lapfit)), fy = nu*sqrt(sample_size)/sqrt(2)*exp(cy)*dnorm(sqrt(sample_size)*(-y_scaled-b)), laplaceVar = 1/sample_size - 4/sample_size^2/nu^2*((Fy+Gy)*fy-2*Fy*Gy)/(Fy+Gy)^2, laplaceVar = if_else(is.nan(laplaceVar),1/sample_size, laplaceVar))
      }
    )
  )


