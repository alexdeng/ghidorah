library(R6)

#' R6 class for a EB method
#' to create an instance, use simpleSplitting$new(trainingSet) where trainingSet is a tibble with required schema to be checked in the initialize function
simpleSplitting = R6Class("RwES-LinearReg",
  public = list(
    includeVar = FALSE,
    initialize= function(trainingSet){
      required = c("y_scaled_A","y_scaled_B","sample_size_A")
      assertthat::assert_that(length(intersect(required,colnames(trainingSet))) == length(required), msg = "Input does not have all the required columns: y_scaled_A, y_scaled_B, sample_size_A!")
      private$training = trainingSet %>% select(required)
    },
    train = function(){
      enriched = (private$training) %>% (private$enrichdata)
      private$model = lm(y_scaled_O ~ y_scaled, data = enriched, weights = enriched$sample_size_A/1000)
      invisible(self)
    },
    predict = function(newdata){
      required = c("y_scaled","scale","sample_size")
      assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
      scaledPred = predict(private$model, newdata = newdata)
      if(self$includeVar){
        return(list(fit =scaledPred * newdata$scale, varfit = 1/newdata$sample_size*newdata$scale^2)) # need to change varfit
      }
      scaledPred * newdata$scale
    },
    printModel = function() print(private$model)
  ),                        
  private = list(
    training = NULL,
    model = NULL,
    enrichdata = function(df){
      df %>% rename(y_scaled=y_scaled_A, y_scaled_O=y_scaled_B)
    }
  ),
  active = list(
    trained = function(){
      !is.null(private$model)
    }
  )
)

