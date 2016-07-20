#' @title Convert set of PNGs to PDF
#' @description Converts PNGs to PDF
#' @param pngs List of names (or with wildcards) of png files
#' @param pdfname Output PDF filename
#' @param extra.opts Options passed to \code{\link{im.convert}}
#' @return Name of PDF output file
#' @export
#' @importFrom animation ani.options im.convert
png_converter <- function(pngs,
                          pdfname,
                          extra.opts = "-density 300"){

  ######################################
  # Turn off the animation options
  ######################################
  aniopts = ani.options()
  ani.options(autobrowse = FALSE)
  ani.options(interval = 0)
  im.convert(pngs,
             output = pdfname,
             extra.opts = extra.opts)
  ######################################
  # Reinstate
  #######################################
  ani.options(aniopts)
  pdfname
}