#' @title Row Centering and Scaling
#' @description This function is a fast way to scale rows of a matrix
#' @param x numeric matrix
#' @param center Should the rows be centered
#' @param scale Should the rows be scaled
#' @param add_attr Add the center/scale attributes
#' @param rows A vector indicating subset of rows to operate over.
#' If \code{NULL}, no subsetting is done.
#' @param cols A vector indicating subset of columns to operate over.
#' If \code{NULL}, no subsetting is done.
#' @param na.rm If \code{TRUE}, missing values are removed first, otherwise not.
#' @param ... Arguments to pass to \code{\link{rowSds}}
#' @return Matrix of centered/scaled values
#' @export
#' @importFrom matrixStats rowSds
#' @examples
#' x = matrix(rnorm(10*1000, mean = 4, sd = 5), ncol = 10)
#' cx = rowScale(x)
#' all(abs(rowMeans(cx)) < 1e-8)
#' all(abs(matrixStats::rowSds(cx) - 1) < 1e-8)
rowScale = function(x,
                    center = TRUE,
                    scale = TRUE,
                    add_attr = TRUE,
                    rows = NULL,
                    cols = NULL,
                    na.rm = TRUE, ...){

  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }

  ################
  # Get the column means
  ################
  cm = rowMeans(x, na.rm = na.rm)
  ################
  # Get the column sd
  ################
  if (scale) {
    csd = rowSds(x, center = cm, na.rm = na.rm, ...)
  } else {
    # just divide by 1 if not
    csd = rep(1, length = length(cm))
  }
  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }
  x =  (x - cm) / csd
  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
  }
  return(x)
}