#' @title Get the AUC and pAUC for a prediction
#' @description Get the area under the curve (AUC) and the partial AUC (pAUC)
#' for a prediction object
#' @param prediction.obj An object of class \code{\link{prediction}}.
#' @param fpr.stop Limit of false positive rate.
#'
#' @return Depending on the number of predictions, either a vector of length 2 or
#' a matrix of 2 columns, with the first column being the AUC, the second - pAUC
#' @export
#' @importFrom ROCR performance
auc_pauc = function(prediction.obj, fpr.stop = 0.01){
  auc = unlist(performance(prediction.obj, "auc")@y.values)
  pauc = unlist(performance(prediction.obj, "auc",
                            fpr.stop = fpr.stop)@y.values) / fpr.stop
  res = c(auc = auc, pauc = pauc)
  if (length(auc) > 1) {
    res = cbind(auc = auc, pauc = pauc)
  }
  return(res)
}

