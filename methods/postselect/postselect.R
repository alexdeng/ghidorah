library(R6)
source('./methods/postselect/selective.R')

postSelectZCut = function(zcut){
  R6Class("postSelect(|z|>${zcut})" %>% stringr::str_interp(list(zcut=zcut)),
          public = list(
            zcut = zcut,
            includeVar = TRUE,
            asymCI = TRUE,
            # post selection inference do not require training. Init by training set just to make the same interface with others
            initialize= function(trainingSet){
              required = c("y", "scale", "sample_size")
              assertthat::assert_that(length(intersect(required,colnames(trainingSet))) == length(required), msg = "Input does not have all the required columns: y, scale, sample_size!")
              private$training = trainingSet %>% select(required)
            },
            # post selection inference do not require training. Init by training set just to make the same interface with others
            train = function(){
              invisible(self)
            },
            predict = function(newdata){
              required = c("y","scale","sample_size")
              assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y,  scale, sample_size!")
              n = length(newdata$y)
              preds = sapply(1:n, function(i){selectiveMLE(newdata$y[i], newdata$scale[i]/sqrt(newdata$sample_size[i]),"both", zcut=self$zcut)})
              if(self$includeVar){
                cis = vapply(1:n, function(i){selectiveCI(newdata$y[i], newdata$scale[i]/sqrt(newdata$sample_size[i]),"both", zcut=self$zcut, cisize = 0.95)}, FUN.VALUE = c(0,0))
                ci = list(lower = cis[1,],upper=cis[2,])
                #implied variance (not symmetric, just to compare CI width)
                varfit = ((ci$upper-ci$lower)/1.96)^2
                return(list(fit = preds, varfit = varfit, ciupper = ci$upper, cilower = ci$lower))
                #list(fit = preds, varfit = varfit, FQfit = newdata$y_scaled*newdata$scale, pMove = pMove, pval = pval)
              }
              return(preds)
            },
            printModel = function() print("Post selection inference with conditinal MLE.")
          ),                        
          private = list(
            training = NULL,
            model = NULL
          ),
          active = list(
            trained = function(){
              TRUE
            }
          )
  )
}