library(R6)
# Baseline unajusted method: use y to predict its mean without any shrinkage, use standard CI 


baselineUnadjust = R6Class("baselineUnadjust",
   public = list(
     includeVar = FALSE,
     # baseline inference do not require training. Init by training set just to make the same interface with others
     initialize= function(trainingSet){
       required = c("y","scale", "sample_size")
       assertthat::assert_that(length(intersect(required,colnames(trainingSet))) == length(required), msg = "Input does not have all the required columns: y, scale, sample_size!")
       private$training = trainingSet %>% select(required)
     },
     train = function(){
       invisible(self)
     },
     predict = function(newdata){
       required = c("y","scale","sample_size")
       assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y, scale, sample_size!")
       if(self$includeVar){
         return(list(fit = newdata$y, varfit = with(newdata, scale^2/sample_size)))
       }
       return(newdata$y)
     }
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



