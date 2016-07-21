#' @title Make Data.frame of images
#' @description Take the filenames and make an analytic data.frame
#' @param fnames Filenames to put into data.frame
#' @param mask Mask to subset the images
#' @param y_fname Name of gold standard
#'
#' @return List of mask and data.frame
#' @importFrom stats model.matrix
#' @export
make_df = function(fnames, mask, y_fname = NULL) {
    mask = check_nifti(mask)

    dd = llply(fnames, function(x){
      img = readnii(x)
      vals = img[ mask == 1 ]
      nonfinite = !is.finite(vals)
      if (any(nonfinite)) {
        vals[nonfinite] = 0
      }
      vals
    }, .progress = "text")


    df = as.data.frame(dd)
    if (!is.null(y_fname)) {
      y = readnii(y_fname)
      df$Y = y[ mask == 1 ]
    }

    # x_names = colnames(df)

    if ("gclass" %in% colnames(df)) {
      df$gclass = factor(df$gclass)
    }
    if ("class" %in% colnames(df)) {
      df$class = factor(df$class)
    }
    df$ants_seg = factor(df$ants_seg)
    df$flair_ants_seg = factor(df$flair_ants_seg)

    cn = c("class",
           "gclass",
           "ants_seg")

    mm = model.matrix(~ . - 1,
                      data = df[, cn])

    df = cbind(df, mm)

    df$class = NULL
    df$gclass = NULL
    df$ants_seg = NULL


    # L = list(mask = mask,
    #      df = df)
    # return(L)
    return(df)
}