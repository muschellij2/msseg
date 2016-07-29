#' @title Convert set of PNGs to PDF
#' @description Converts PNGs to PDF
#' @param pngs List of names (or with wildcards) of png files
#' @param pdfname Output PDF filename
#' @param extra.opts Options passed to \code{\link{im.convert}}
#' @param interval interval of delay, passed to \code{\link{ani.options}}
#' @param autobrowse passed to \code{\link{ani.options}}
#' @return Name of PDF output file
#' @export
#' @importFrom animation ani.options im.convert
png_converter <- function(pngs,
                          pdfname,
                          extra.opts = "-density 300",
                          interval = 0,
                          autobrowse = FALSE){

  ######################################
  # Turn off the animation options
  ######################################
  aniopts = ani.options()
  ani.options(autobrowse = autobrowse)
  ani.options(interval = interval)
  im.convert(pngs,
             output = pdfname,
             extra.opts = extra.opts)
  ######################################
  # Reinstate
  #######################################
  ani.options(aniopts)
  pdfname
}