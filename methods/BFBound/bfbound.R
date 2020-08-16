library(R6)

pvalBFBound = function(pval){
  bfb = 1/(-exp(1)*pval*log(pval))
  bfb[pval>exp(-1)] = exp(-1)/pval[pval>exp(-1)] #e^-1/p is better than 1 . These are noninteresting case anyway
  bfb
}

localH1Bound = function(...){
  vars = rlang::dots_list(...)
  if(length(vars)>1){
    ph1 = vars[[1]]
    ph0 = vars[[2]]
    priorOdds = ph1/ph0
    nametag = paste0("LocalH1(",ph1,":",ph0,")")
  } else{
    priorOdds = vars[[1]]
    nametag = paste0("LocalH1(pOdds=",priorOdds,")")
  }
  R6Class(nametag, 
    public = list(
      priorOdds = priorOdds,
      includeVar = TRUE,
      initialize= function(trainingSet) {}, # no action needed
      train = function(){invisible(self)}, 
      predict = function(newdata){
        required = c("y_scaled","scale","sample_size")
        assertthat::assert_that(length(intersect(required,colnames(newdata))) == length(required), msg = "Input does not have all the required columns: y_scaled, scale, sample_size!")
        tstat = abs(newdata$y_scaled*sqrt(newdata$sample_size))
        pval = 2*pnorm(tstat, 0, 1, lower=FALSE)
        BF1v0 = pvalBFBound(pval)
        BF1v0[is.nan(BF1v0)]=1e6
        postOdds = BF1v0*self$priorOdds
        pMove = postOdds/(1+postOdds)
        #print(pMove)
        preds = newdata$y_scaled*pMove*newdata$scale
        if(self$includeVar){
          varfit = (pMove*1/newdata$sample_size + pMove*newdata$y_scaled^2 - pMove^2*newdata$y_scaled^2)*newdata$scale^2
          fqvar = 1/newdata$sample_size*newdata$scale^2
          ind = varfit>fqvar #faster than ifelse
          varfit[ind] = fqvar[ind] # cap at FQ variance.
          return(list(fit = preds, varfit = varfit, FQfit = newdata$y_scaled*newdata$scale, pMove = pMove, pval = pval))
        }
        return(list(fit = preds, FQfit = newdata$y_scaled*newdata$scale, pMove = pMove, pval = pval))      
      }
    )
  )
}




