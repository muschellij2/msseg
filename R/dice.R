#' @title Calculate Dice Coefficient
#' @description Calculate Dice Coefficient/Similarity Index from a prediction object
#' @param prediction.obj object of class \code{\link{prediction-class}} from
#' \code{\link{prediction}}
#' @return Object of class \code{\link{performance}}
#' @export
#' @importFrom methods new
dice <- function(prediction.obj){
  if (class(prediction.obj) != "prediction") {
    stop(paste(
      "Wrong argument types: First argument",
      "must be of type",
      "'prediction'"))
  }
  argnames <- c()
  x.values <- list()
  y.values <- list()
  for (i in 1:length(prediction.obj@predictions)) {
    fp = prediction.obj@fp[[i]]
    tp = prediction.obj@tp[[i]]
    fn = prediction.obj@fn[[i]]
    tn = prediction.obj@tn[[i]]
    cutoffs = prediction.obj@cutoffs[[i]]
    meas_dice = 2 * tp / (2*tp + fp + fn)
    x.values <- c(x.values, list(cutoffs))
    y.values <- c(y.values, list(meas_dice))
  }
  if (!(length(x.values) == 0 ||
        length(x.values) == length(y.values))) {
    stop("Consistency error.")
  }
  return(new("performance", x.name = "cutoff",
             y.name = "dice",
             alpha.name = "none",
             x.values = x.values, y.values = y.values,
             alpha.values = list())
  )
}