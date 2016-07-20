
#' @title Trimmed Mean and Standard Deviation
#' @description Calculates both trimmed mean and standard deviation from a vector
#' @param x numeric vector
#' @param trim value between 0 and 0.5, percent trim from each side of the distribution.
#' Calculates quantiles for \code{trim} and \code{1-trim} and sets values outside to NA
#'
#' @return Vector of length 2, with mean/sd named elements
#' @export
#' @examples
#' x = rnorm(100)
#' trim_mean_sd(x)
#' @importFrom stats quantile sd
trim_mean_sd <- function(x, trim = 0.2){
  qtrim <- quantile(x,
                    c(trim, 0.5, 1 - trim),
                    na.rm = FALSE)
  xbot <- qtrim[1]
  xtop <- qtrim[3]
  x[x < xbot] <- NA
  x[x > xtop] <- NA

  mn = mean(x, na.rm = TRUE)
  s = sd(x, na.rm = TRUE)
  c(mean = mn, sd = s)
}