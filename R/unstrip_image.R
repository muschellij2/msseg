#' @title Replace Dropped Dimensions from Stripping
#' @description Takes in a _strip segmentation and makes it same size as FLAIR
#' @param strip_image image with dropped dimensions
#' @param outdir Location of \code{FLAIR_MALF_Brain_Mask.nii.gz} and
#' \code{FLAIR_noneck.nii.gz}
#' @return Nifti object of same dimensions as FLAIR
#' @export
#' @importFrom neurobase replace_dropped_dimensions
#' @importFrom extrantsr check_ants
unstrip_image = function(strip_image,
         outdir = "."){

  ##############################################
  # First Quick Pass
  ##############################################
  img = check_nifti(strip_image)

  # mask_fname = file.path(outdir,
  #                        "FLAIR_MALF_Brain_Mask_strip.nii.gz")
  # stopifnot(file.exists(mask_fname))
  # mask = readnii(mask_fname)

  orig_mask = file.path(outdir,
                        "FLAIR_MALF_Brain_Mask.nii.gz")
  stopifnot(file.exists(orig_mask))
  mask_full = readnii(orig_mask)

  dd_all = dropEmptyImageDimensions(
    mask_full,
    keep_ind = TRUE)

  img = replace_dropped_dimensions(img = img,
                                   inds = dd_all$inds,
                                   orig.dim = dim(mask_full))

  ##############################################
  # First Quick Pass
  ##############################################
  nn = file.path(outdir, "FLAIR_noneck.nii.gz")
  stopifnot(file.exists(nn))
  x = check_ants(nn)
  nn = readnii(nn)

  gm = getMask(x, cleanup = 0)
  rm(list = c("x")); gc();
  x = ants2oro(gm)
  rm(list = c("gm")); gc();
  dd = dropEmptyImageDimensions(x, keep_ind = TRUE)

  img = replace_dropped_dimensions(img = img,
                                   inds = dd$inds,
                                   orig.dim = dim(nn))
  return(img)
}