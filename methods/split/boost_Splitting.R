library(R6)
library(xgboost)
#gradient boosting tree
bstSplitting = R6Class("bstSplitting",
  public = list(
   initialize= function(trainingSet){
      required = c("y_scaled_A","y_scaled_B","sample_size_A")
      assertthat::assert_that(length(intersect(required,colnames(trainingSet))) == length(required), msg = "Input does not have all the required columns: y_scaled_A, y_scaled_B, sample_size_A!")
     private$training = trainingSet %>% select(required)
   },
   train = function(verbose=1){
     enriched = (private$training) %>% (private$enrichTrainingData)
     private$model <- xgboost::xgb.train(data = enriched, booster = "gbtree", max.depth = 2, eta = 0.6, nthread = 4, nrounds = 6, objective = "reg:squarederror",verbose = verbose)
     invisible(self)
   },
   predict = function(newdata){
     required = c("y_scaled","scale","sample_size")
     assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
     scaledPred = predict(private$model, newdata = newdata %>% (private$enrichTestingData))
     scaledPred * newdata$scale
   },
   printModel = function() print(private$model),
   getModel = function() (private$model)
  ),                        
  private = list(
   training = NULL,
   model = NULL,
   enrichTrainingData = function(df){
     tmp = df %>% rename(y_scaled=y_scaled_A, y_scaled_O=y_scaled_B) %>% mutate(sample_size = sample_size_A)  # sample size for split is only a half
     y = tmp$y_scaled_O
     x = tmp %>% (private$transformFeature) %>%  select(y_scaled, sample_size) %>% data.matrix
     xgboost::xgb.DMatrix(data = x, label = y, weight = tmp$sample_size/1000)
   },
   enrichTestingData = function(df){
     df %>% (private$transformFeature) %>% select(y_scaled, sample_size) %>% data.matrix
   },
   transformFeature = function(df) df  # more teature later
  ),
  active = list(
   trained = function(){
     !is.null(private$model)
   }
  )
)

#linear gradient boosting
lbsSplitting = R6Class("lbsSplitting",
  public = list(
   initialize= function(trainingSet){
      required = c("y_scaled_A","y_scaled_B","sample_size_A")
      assertthat::assert_that(length(intersect(required,colnames(trainingSet))) == length(required), msg = "Input does not have all the required columns: y_scaled_A, y_scaled_B,, sample_size_A!")
     private$training = trainingSet %>% select(required)
   },
   train = function(verbose=0){
     enriched = (private$training) %>% (private$enrichTrainingData)
     private$model <- xgboost::xgb.train(data = enriched, booster = "gblinear", nthread = 4, nrounds = 6, objective = "reg:squarederror",verbose = verbose)
     invisible(self)
   },
   predict = function(newdata){
     required = c("y_scaled","scale","sample_size")
     assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
     scaledPred = predict(private$model, newdata = newdata %>% (private$enrichTestingData))
     scaledPred * newdata$scale
   },
   printModel = function() print(private$model),
   getModel = function() (private$model)
  ),                        
  private = list(
   training = NULL,
   model = NULL,
   enrichTrainingData = function(df){
     tmp = df %>% rename(y_scaled=y_scaled_A, y_scaled_O=y_scaled_B) %>% mutate(sample_size = sample_size_A)  # sample size for a split is a half
     y = tmp$y_scaled_O
     x = tmp %>% (private$transformFeature) %>%  select(y_scaled, sample_size) %>% data.matrix
     xgboost::xgb.DMatrix(data = x, label = y, weight = tmp$sample_size/1000)
   },
   enrichTestingData = function(df){
     df %>% (private$transformFeature) %>% select(y_scaled, sample_size) %>% data.matrix
   },
   transformFeature = function(df) df  # more teature later
  ),
  active = list(
   trained = function(){
     !is.null(private$model)
   }
  )
)
