#' @title Predict environmental effects
#'
#' @description A function to train a random forest model to predict environmental effects defined from a MET mixed model analysis with ECs as predictors.
#'
#' @param train.ECs A data frame of weather and/or soil ECs for observed environments as output from the [get.W.ECs()] or [get.S.ECs()] functions.
#' @param new.ECs A data frame of weather and/or soil ECs for new environments as output from the [get.W.ECs()] or [get.S.ECs()] functions.
#' @param E.effs A data frame of several or a vector of a single environmental effects parameters fitted from a multi-environmnet trial analysis mixed model.
#' Latent environmental effect factor loadings that decompose GxE can be defined as described by [Smith et al (2021)](https://doi.org/10.3389/fpls.2021.737462)
#' @param verbose Logical. Should progress be printed? Default = TRUE.
#' @param ntree Number of decision trees to use in random forest models. Default = 1000 trees.
#' @param ncores Number (integer) of cores to use for parallel processing. The default (`NULL`) will use the one less than the maximum available.
#'
#' @details
#' Random forest models are fitted with one 3rd of feature variables samples at each split and a minimum end node size of 5.
#' 
#' @returns A data frame of environmental effect predictions for the new environments with environments as rows and environmental effect variates as columns.
#'
#' @references
#' Smith, A., Norman, A., Kuchel, H., & Cullis, B. (2021). Plant Variety Selection Using Interaction Classes Derived From Factor Analytic Linear Mixed Models:
#'      Models With Independent Variety Effects. Frontiers in Plant Science, 12. <https://doi.org/10.3389/fpls.2021.737462>
#'
#' @author Nick Fradgley
#'
#' @export

pred.env.effs <- function(train.ECs, new.ECs, E.effs, verbose=TRUE, ntree=1000, ncores = NULL) {
  
  E.effs <- as.data.frame(E.effs)
  if (sum(!rownames(train.ECs) == rownames(E.effs)) > 0) {
    warning("Row names for ECs and E.effs do not match\n")
  }

  if (sum(!colnames(train.ECs) %in% colnames(new.ECs)) > 0) {
    warning("Variables in train.ECs but not in new.ECs\n")
  }
  
  if(sum(!apply(E.effs,2,is.numeric))>0){
    warning("E.effs non-numeric\n")
  }
  
  if (is.null(ncores)) {
    ncores <- parallel::detectCores()-1
  }
  
  if (verbose) cat("Training RF models\n")
  rfmods <- apply(E.effs, 2, function(y) {
    if (verbose) cat("|")
    ranger::ranger(x = train.ECs, y = y, num.trees = ntree,mtry = round(ncol(train.ECs)/3),min.node.size =  5,num.threads = ncores )
    
  })

  if (verbose) cat("\nPredicting to new envs\n")
  preds <- sapply(rfmods, function(x) {
    if (verbose) cat("|")
    predict(x, data = new.ECs[, colnames(train.ECs)])$predictions
  })
  preds <- as.data.frame(preds)
  rownames(preds)<-rownames(new.ECs)
  return(preds)
}


