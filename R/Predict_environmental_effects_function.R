#' @title Predict environmental effects
#'
#' @description A function to train a random forest model to predict environmental effects defined from a MET mixed model analysis with ECs as predictors.
#'
#' @param train.ECs A data frame of weather and/or soil ECs for observed environments as output from the [get.W.ECs()] or [get.S.ECs()] functions.
#' @param new.ECs A data frame of weather and/or soil ECs for new environments as output from the [get.W.ECs()] or [get.S.ECs()] functions.
#' @param E.effs A data frame of several or a vector of a single environmental effects parameters fitted from a multi-environmnet trial analysis mixed model.
#' Latent environmental effect factor loadings that decompose GxE can be defined as described by [Smith et al (2021)](https://doi.org/10.3389/fpls.2021.737462)
#' @param verbose Logical. Should progress be printed? Default = TRUE.
#' @param ntree Number of decision trees to use in random forest models.
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

pred.env.effs <- function(train.ECs, new.ECs, E.effs, verbose=TRUE, ntree=500) {
  
  E.effs <- as.data.frame(E.effs)
  if (sum(!rownames(train.ECs) == rownames(E.effs)) > 0) {
    print("Row names for ECs and E.effs do not match\n")
  }

  if (sum(!colnames(train.ECs) %in% colnames(new.ECs)) > 0) {
    print("Variables in train.ECs but not in new.ECs\n")
  }
  cat("Training RF models\n")
  rfmods <- apply(E.effs, 2, function(x) {
      cat("|")
    randomForest::randomForest(x = train.ECs, y = x, ntree = ntree)
  })

  cat("Training RF models\n")
  preds <- sapply(rfmods, function(x) {
    predict(x, newdata = new.ECs[, colnames(train.ECs)])
  })
  preds <- as.data.frame(preds)
  return(preds)
}
