
#' @title Get the Maximum Dice Coefficient and Cutoff
#' @description Get the Dice coefficient and the cutoff for that max Dice
#' from a prediction object
#' @param prediction.obj An object of class \code{\link{prediction}}.
#'
#' @return Depending on the number of predictions, either a vector of length 2 or
#' a matrix of 2 columns, with the first column being the dice, the second - cutoff
#' @export
get_max_dice = function(prediction.obj){
  mydice = dice(prediction.obj)
  # cutoffs = mydice@x.values[[1]]
  # mydice = mydice@y.values[[1]]
  res = mapply(function(cutoffs, dices){
    names(cutoffs) = names(dices) = NULL
    max_dice_ind = which.max(dices)
    max_dice = dices[max_dice_ind]
    max_dice_cutoff = cutoffs[max_dice_ind]
    c(max_dice = max_dice,
      max_dice_cutoff = max_dice_cutoff)
  }, mydice@x.values, mydice@y.values, SIMPLIFY = TRUE)
  res
}