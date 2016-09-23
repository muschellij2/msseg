
#' @title Trimmed Z-score Image
#' @description Creates a normalized image using z-scores
#' @param img Object of class \code{nifti}
#' @param mask Binary object of class \code{nifti}.  If missing, whole image used.
#' @param remask Should the result be masked after z-scoring (all outside set to 0)
#' @param trim The trim passed to \code{\link{trim_mean_sd}}
#' @param center Logical indicating whether the value should be centered by
#' trimmed mean
#' @param scale Logical indicating whether the value should be scaled by
#' trimmed sd
#' @return Object of class \code{nifti}
#' @export
#' @importFrom neurobase niftiarr check_nifti check_mask_fail mask_img
#' @examples
#' nim = oro.nifti::nifti(img = array(rnorm(1000), dim = c(10, 10, 10)))
#' mask = nim > 0
#' trimmed_z(nim, mask)
trimmed_z <- function(img, mask = NULL,
                  remask = FALSE, trim = 0.2,
                  center = TRUE,
                  scale = TRUE){
  img = check_nifti(img)
  if (is.null(mask)) {
    mask = niftiarr(img, 1)
  }
  check_mask_fail(mask, allow.NA = FALSE)

  x = img[ mask == 1 ]
  tvals = trim_mean_sd(x, trim = trim)
  mn = tvals['mean']
  s = tvals['sd']
  img = (img-mn)/s
  if (remask) {
    img = mask_img(img, mask)
  }
  return(img)
}


#' @title Trimmed Ratio Image
#' @description Creates a normalized image using the image divided by a
#' trimmed standard deviation, but not centered
#' @param ... Arguments passed to \code{\link{trimmed_z}}
#' @return Object of class \code{nifti}
#' @export
#' @examples
#' nim = oro.nifti::nifti(img = array(rnorm(1000), dim = c(10, 10, 10)))
#' mask = nim > 0
#' trimmed_z(nim, mask)
trimmed_ratio = function(...){
  img = trimmed_z(..., center = FALSE)
  return(img)
}