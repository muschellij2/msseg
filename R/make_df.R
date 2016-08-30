#' @title Make Data.frame of images
#' @description Take the filenames and make an analytic data.frame
#' @param fnames Filenames to put into data.frame
#' @param mask Mask to subset the images
#' @param y_fname Name of gold standard
#' @param verbose print diagnostic messages
#'
#' @return List of mask and data.frame
#' @importFrom stats model.matrix
#' @importFrom dplyr bind_cols
#' @export
make_df = function(fnames, mask,
                   y_fname = NULL,
                   verbose = TRUE) {
    mask = check_nifti(mask)

    if (verbose) {
      message("reading in data")
    }
    df = llply(fnames, function(x){
      img = readnii(x)
      vals = img[ mask == 1 ]
      nonfinite = !is.finite(vals)
      if (any(nonfinite)) {
        vals[nonfinite] = 0
      }
      vals
    }, .progress = "text")

    if (verbose) {
      message("Converting to data.frame")
    }
    df = as.data.frame(df)
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
    if (verbose) {
      message("Making model.matrix output")
    }
    mm = model.matrix(~ . - 1,
                      data = df[, cn])
    mm = as.data.frame(mm)
    if (verbose) {
      message("Making model.matrix output")
    }
    df$class = NULL
    df$gclass = NULL
    df$ants_seg = NULL

    if (verbose) {
      message("Binding Columns")
    }
    for (icol in colnames(mm)) {
      df[, icol] = mm[, icol]
    }
    # df = bind_cols(df, mm)

    # L = list(mask = mask,
    #      df = df)
    # return(L)
    return(df)
}