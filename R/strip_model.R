#' @title Strip Model to Reduce Memory Footprint
#' @description Strips elements from a \code{glm} or \code{lm} object to reduce size.
#' Allows for printing/showing, but \code{summary} must be run before stripping
#' @param model Object of class \code{glm} or \code{lm}
#' @return Object of class \code{glm} or \code{lm}, but with some elements removed
#' @export
#' @examples
#' ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
#' trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
#' group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
#' weight <- c(ctl, trt)
#' lm.D9 <- lm(weight ~ group)
#' object.size(lm.D9)
#' object.size(strip_model(lm.D9))
strip_model = function(model){
  model$y = c()
  model$model = c()
  model$residuals = c()
  model$fitted.values = c()
  model$effects = c()
  model$qr$qr = c()
  model$linear.predictors = c()
  model$weights = c()
  model$prior.weights = c()
  model$data = c()
  attr(model$terms,".Environment") = c()
  attr(model$formula,".Environment") = c()
  model
}